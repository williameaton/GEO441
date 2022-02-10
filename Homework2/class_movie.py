import struct


with open('slices.out', mode='rb') as file: # b is important -> binary
    slices = file.read()

a = struct.unpack("iiiii", slices[:20])

print()