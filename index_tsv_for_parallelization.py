#%%
import struct, math, re, sys

_bgzf_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
_bgzf_header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"
_bgzf_magic = b"\x1f\x8b\x08\x04"


def index_txt_file(path_file, path_save = None, comment_indicator = None):
    """_summary_
    Generate Text Index (TXI) file for Indexing line with line break (\n).
    Format:
        n_comment_line (uint32)
        n_data_line (uint32)
        offset_header (uint32)
        offset_data_line (uint32[n_data_line])

    Args:
        path_file (str): Path to TXT file for indexing
        path_save (str): Path to save TXI file. Defaults to {path_file}.txi
        comment_indicator (str, optional): Indicator (First character) of comment lines to ignore. Defaults to None.
    """
    dict_file_offsets = get_line_info(path_file, comment_indicator)
    if path_save == None:
        path_save = f"{path_file}.txi"
    write_txi(path_save, dict_file_offsets)    

def get_line_info(path_file, comment_indicator):
    file_handler = open(path_file, 'r')
    
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
                current_offset = file_handler.tell()
                break
    dict_file_offsets["Header"].append(current_offset)
    header = file_handler.readline()
    current_offset = file_handler.tell()
    
    while 1:
        line = file_handler.readline()
        if not line:
            break
        dict_file_offsets["Data"].append(current_offset)
        current_offset = file_handler.tell()
    
    file_handler.close()
    return dict_file_offsets    

def write_txi(path_save, dict_file_offsets):
    file_handler = open(path_save, "wb")
    
    n_comment_line = len(dict_file_offsets["Comment"])
    n_data_line = len(dict_file_offsets["Data"])
    
    list_content = [n_comment_line, n_data_line, dict_file_offsets["Header"][0]]
    list_content.extend(dict_file_offsets["Data"])
    list_content_binary = list(map(lambda val: struct.pack("<I", val), list_content))
    
    binary_index_data = b''.join(list_content_binary)
    file_handler.write(binary_index_data)
    file_handler.flush()
    file_handler.close()
    
def read_txt_index(path_index):
    file_handler = open(path_index, 'rb')
    n_comment_line, n_data_line, header_offset = struct.unpack("<III", file_handler.read(12))
    
    data_line_offsets = file_handler.read()
    list_line_offset = [struct.unpack("<I", data_line_offsets[ind*4:ind*4+4])[0] for ind in range(n_data_line)]
    
    dict_file_offsets = {
        "N_line":{"Comment":n_comment_line, "Data":n_data_line},
        "Header": [header_offset],
        "Data":list_line_offset
    }
    file_handler.close()
    return dict_file_offsets
# %%
