#%%
from Bio import bgzf
import struct, math, re, sys
import numpy as np
from joblib import Parallel, delayed
from itertools import chain

_bgzf_magic = b"\x1f\x8b\x08\x04"
_bgzf_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"


def index_bgzipped_txt_file(path_file, path_save = None, comment_indicator = None, parallel = 1):
    """_summary_
    Generate Block Gzipped Text Index (.gz) file for Indexing line with line break (\n).
    Format:
        n_comment_line (uint32)
        n_data_line (uint32)
        virtual_offset_header (uint32)
        virtual_offset_data_line (uint32[n_data_line])

    Args:
        path_file (_type_): _description_
        path_save (_type_, optional): _description_. Defaults to None.
        comment_indicator (_type_, optional): _description_. Defaults to None.
    """
    # dict_file_offsets = get_line_info(path_file, comment_indicator)
    dict_line_offsets = get_header_line_info(path_file, comment_indicator)
    whole_linestart_virtual_offsets = find_bgzip_data_line_offsets_parallely(path_file, parallel)
    list_data_linestart_virtual_offsets = sorted(set(whole_linestart_virtual_offsets) - set(dict_line_offsets["Comment"]+dict_line_offsets["Header"]))
    
    dict_line_offsets["Data"] = list_data_linestart_virtual_offsets    
    
    if path_save == None:
        path_save = f"{path_file}.gzli"
    write_bzli(path_save, dict_line_offsets)       
    
def get_header_line_info(path_file, comment_indicator):
    file_handler = bgzf.BgzfReader(path_file)
    
    dict_file_offsets = {
        "Comment": [],
        "Header": []
    }
    current_offset = 0
    if comment_indicator:
        while 1:
            line = file_handler.readline()
            if line.startswith(comment_indicator):
                dict_file_offsets['Comment'].append(current_offset)
                current_offset = file_handler.tell()
            else:
                dict_file_offsets["Header"].append(current_offset)
                current_offset = file_handler.tell()
                break
    
    file_handler.close()
    return dict_file_offsets

def write_bzli(path_save, dict_file_offsets):
    file_handler = open(path_save, "wb")
    
    n_comment_line = len(dict_file_offsets["Comment"])
    n_data_line = len(dict_file_offsets["Data"])
    
    list_content_binary = [
        struct.pack("<I", n_comment_line),
        struct.pack("<I", n_data_line),
        struct.pack("<Q", dict_file_offsets["Header"][0])
    ]
    for val in dict_file_offsets["Data"]:
        bsize = val >> 16
        within_block = val ^ (bsize << 16)
        list_content_binary.append(struct.pack("<Q", bsize))
        list_content_binary.append(struct.pack("<H", within_block))
    binary_index_data = b''.join(list_content_binary)
    file_handler.write(binary_index_data)
    file_handler.flush()
    file_handler.close()
    
def read_bgzip_text_index(path_index):
    file_handler = open(path_index, 'rb')
    n_comment_line, n_data_line, header_offset = struct.unpack("<IIQ", file_handler.read(16))
    
    data_line_offsets = file_handler.read()
    list_line_offset = [struct.unpack("<QH", data_line_offsets[ind*10:ind*10+10]) for ind in range(n_data_line)]
    
    dict_file_offsets = {
        "N_line":{"Comment":n_comment_line, "Data":n_data_line},
        "Header": [header_offset],
        "Data":list_line_offset
    }
    file_handler.close()
    return dict_file_offsets

def find_bgzip_data_line_offsets_parallely(path_file, n_parallel):
    bytes_of_bgzip_eof = get_end_of_file_position(path_file)
    list_search_start_pos = list(map(int, np.linspace(0, bytes_of_bgzip_eof, n_parallel+1)[:-1]))
    
    file_handler = open(path_file, "rb")
    list_block_start_offsets = list(map(lambda bytes_start_search: search_nearest_bgzip_block(file_handler, bytes_start_search), list_search_start_pos))
    file_handler.close()
    list_block_end_offsets = list_block_start_offsets[1:]+[bytes_of_bgzip_eof]
    
    with Parallel(n_parallel) as parallel:
        list_list_lineend_virtual_offsets = parallel(delayed(get_line_end_virtual_offsets_between_two_offsets)(
            path_file,
            bytes_block_start,
            bytes_block_end
        )for (bytes_block_start, bytes_block_end) in zip(list_block_start_offsets, list_block_end_offsets))
    list_whole_line_end_virtual_offsets = list(chain(*list_list_lineend_virtual_offsets))
        
    sorted_whole_lineend_virtual_offsets = sorted(set(list_whole_line_end_virtual_offsets))
    whole_linestart_virtual_offsets = [0] + sorted_whole_lineend_virtual_offsets[:-1]
    return whole_linestart_virtual_offsets

def get_end_of_file_position(path_file):
    file_object = open(path_file, "rb")
    filesize = file_object.seek(0, 2)
    filesize_except_eof = filesize - len(_bgzf_eof)
    file_object.seek(filesize_except_eof)
    eof = file_object.read(len(_bgzf_eof))
    file_object.close()
    assert eof == _bgzf_eof, "Block gzip file does not ends with 'end of file' context"
    return filesize_except_eof

def search_nearest_bgzip_block(handler, coffset_start_search, wsize_search = 5000, step_search = 4900):
    handler.seek(coffset_start_search)
    curr_coffset = coffset_start_search
    
    while 1:
        handler.seek(curr_coffset)
        dat_search = handler.read(wsize_search)
        try:
            ind_bgzf_magic = dat_search.index(_bgzf_magic)
            curr_coffset += ind_bgzf_magic
            handler.seek(curr_coffset)
            dat_search = handler.read(wsize_search)
            if validate_bgzip_block_header(dat_search):
                break
            else:
                curr_coffset += 4
                continue
        except Exception:
            curr_coffset += step_search
            continue
    return curr_coffset

def validate_bgzip_block_header(data):
    check_gid1 = data[0] == 31
    check_gid2 = data[1] == 139
    check_gcm = data[2] == 8
    check_gflag = data[3] == 4
    check_si1 = data[12] == 66
    check_si2 = data[13] == 67
    check_slen = struct.unpack("<H", data[14:16])[0] == 2
    
    return check_gid1 + check_gid2 + check_gcm + check_gflag + check_si1 + check_si2 + check_slen == 7 

def get_line_end_virtual_offsets_between_two_offsets(path_file, bytes_search_begin, bytes_search_end):
    virtual_offset_begin = bytes_search_begin << 16
    virtual_offset_end = bytes_search_end << 16
    
    file_handler = bgzf.BgzfReader(path_file)
    file_handler.seek(virtual_offset_begin)
    
    list_lineend_offsets = list()
    while 1:
        _ = file_handler.readline()
        end_offset = file_handler.tell()
        list_lineend_offsets.append(end_offset)
        if end_offset >= virtual_offset_end:
            break
    file_handler.close()
    return list_lineend_offsets
# %%
# path_file = "/BiO/Research/Korea10KGenome/Workspace/Yoonsung/ToyData/Jointcalled_vcf.Toyset/korea10K.jointcall.merge156.removeDup.maskNullDepth.remDupSample.chr20.top10000.test.Bio_bgzip.vcf.gz"
# test = get_line_info("/BiO/Research/Korea10KGenome/Workspace/Yoonsung/ToyData/Jointcalled_vcf.Toyset/korea10K.jointcall.merge156.removeDup.maskNullDepth.remDupSample.chr20.top10000.test.Bio_bgzip.vcf.gz", "##")
for chrname in list(map(lambda x: f"chr{x}", range(1, 23))):
    index_bgzipped_txt_file(
        f"/BiO/Research/Korea10KGenome/Results/Korea10K_VCF/Jointcalled_vcf.hard_mask_NA.rem_dup_sample.VQSR.PASS/korea10K.jointcall.merge156.removeDup.maskNullDepth.remDupSample.applyVQSR.PASS.{chrname}.vcf.gz",
        f"/BiO/Access/yoonsung/Research/Custom_BGzip_for_TSV/test.{chrname}.core10.gzli",
        "##",
        10
    )
# %%
