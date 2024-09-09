#%%
from Bio import bgzf
import struct, math, re, sys
import numpy as np

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
    
    dict_file_offsets = get_line_info(path_file, comment_indicator)
    # dict_header_offsets = get_header_line_info(path_file, comment_indicator)
    
    if path_save == None:
        path_save = f"{path_file}.bzli"
    write_bzli(path_save, dict_file_offsets)       
    
def get_line_info(path_file, comment_indicator):
    eof_bytes = get_end_of_file_position(path_file)
    eof_virtual_offset = eof_bytes << 16
    
    file_handler = bgzf.BgzfReader(path_file)
    
    dict_file_offsets = {
        "Comment": [],
        "Header": [],
        "Data": []
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
    
    while 1:
        dict_file_offsets["Data"].append(current_offset)
        line = file_handler.readline()
        current_offset = file_handler.tell()
        if current_offset >= eof_virtual_offset:
            break
        
    file_handler.close()
    return dict_file_offsets

def get_end_of_file_position(path_file):
    file_object = open(path_file, "rb")
    filesize = file_object.seek(0, 2)
    filesize_except_eof = filesize - len(_bgzf_eof)
    file_object.seek(filesize_except_eof)
    eof = file_object.read(len(_bgzf_eof))
    file_object.close()
    assert eof == _bgzf_eof, "Block gzip file does not ends with 'end of file' context"
    return filesize_except_eof

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


# %%
# test = get_line_info("/BiO/Research/Korea10KGenome/Workspace/Yoonsung/ToyData/Jointcalled_vcf.Toyset/korea10K.jointcall.merge156.removeDup.maskNullDepth.remDupSample.chr20.top10000.test.Bio_bgzip.vcf.gz", "##")
# index_bgzipped_txt_file(
#     "/BiO/Research/Korea10KGenome/Workspace/Yoonsung/ToyData/Jointcalled_vcf.Toyset/korea10K.jointcall.merge156.removeDup.maskNullDepth.remDupSample.chr20.top10000.test.Bio_bgzip.vcf.gz",
#     "/BiO/Access/yoonsung/Research/Custom_BGzip_for_TSV/test.bzli",
#     "##"
# )
index_bgzipped_txt_file(
    "/BiO/Research/Korea10KGenome/Results/Korea10K_VCF/Jointcalled_vcf.hard_mask_NA.rem_dup_sample.VQSR.PASS/korea10K.jointcall.merge156.removeDup.maskNullDepth.remDupSample.applyVQSR.PASS.chr1.vcf.gz",
    "/BiO/Access/yoonsung/Research/Custom_BGzip_for_TSV/test.chr1.bzli",
    "##"   
)
# %%
