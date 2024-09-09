#%%
import struct, math, re, sys, zlib

_bgzf_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
_bgzf_header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"
_bgzf_magic = b"\x1f\x8b\x08\x04"




def get_compressed_block_of_data(data, compresslevel = 6):
    if len(data) > 65536:
        raise ValueError(f"{len(data)} Block length > 65536")
    compressor = zlib.compressobj(
        compresslevel, zlib.DEFLATED, -15, zlib.DEF_MEM_LEVEL, 0
    )
    compressed_data = compressor.compress(data) + compressor.flush()
    del compressor
    
    bsize = struct.pack("<H", len(compressed_data) + 25)  # includes -1
    crc = struct.pack("<I", zlib.crc32(data) & 0xFFFFFFFF)
    uncompressed_length = struct.pack("<I", len(data))
    
    compressed_block = _bgzf_header + bsize + compressed_data + crc + uncompressed_length
    return compressed_block

def write_block(file_handler, data):
    if isinstance(file_handler, str):
        file_handler = open(file_handler, "ab")
    compressed_data = get_compressed_block_of_data(data)
    file_handler.write(compressed_data)
    
def write_eof(file_handler):
    if isinstance(file_handler, str):
        file_handler = open(file_handler, "ab")
    file_handler.write(_bgzf_eof)
    file_handler.flush()
    file_handler.close()