import math

def binary2flags(bitstring, base=16):
    
    """
    Convert bitstring to flags in schema with given number of flags; flags are
    labeled A, B, C, D, ...
    """
    
    bitstring = int(bitstring, 2)
    encoding = ''
    while bitstring:
        encoding += chr(ord('A') + (bitstring % base))
        bitstring //= base
    return encoding[::-1]

def flags2binary(flags, base=16):
    
    """
    Convert sequence of flags into corresponding bitstring
    """
    
    return bin(sum(
        (ord(flag) - ord('A')) * (base ** index)
        for index, flag in enumerate(flags[::-1])
    ))[2:]

def getFragmentLength(flag_count, polymer_length):
    
    return math.floor(math.log(flag_count ** polymer_length, 2))

class ChemCoder:
    
    # fixed number of bits to encode the length of the bitstring
    filesize_length = 32

    def __init__(self, flag_count, polymer_length):
        
        self.flag_count = flag_count
        self.polymer_length = polymer_length
    
    def getFragmentLength(self):
        
        return math.floor(math.log(self.flag_count ** self.polymer_length, 2))

    def encode(self, bitstring):
        
        # determine length of the original bitstring
        bitstring_length = len(bitstring)
    
        # prefix bitstring with bitstring length
        bitstring = bin(bitstring_length)[2:].zfill(self.filesize_length) + bitstring
        bitstring_length += self.filesize_length
        
        # compute number of bits per fragment
        fragment_length = self.getFragmentLength()
        fragment_count = math.ceil(bitstring_length / fragment_length)
        
        # determine lower limit of position encoding bit length
        position_length = math.ceil(math.log(fragment_count, 2))
        data_length = fragment_length - position_length
        assert data_length > 0, 'impossible encoding'            
        fragment_count = math.ceil(bitstring_length / data_length)
        
        # update position length until it can encode the combination of position
        # and data fragments
        while 2 ** position_length < fragment_count:
            position_length += 1
            data_length -= 1
            assert data_length > 0, 'impossible encoding'            
            fragment_count = math.ceil(bitstring_length / data_length)
        
        # right pad bitstring to give last fragment the fixed fragment length
        bitstring += '0' * (fragment_count * fragment_length - len(bitstring))
           
        # encode/decode fragments
        fragments = []
        for fragment in range(fragment_count):
            position = bin(fragment)[2:].zfill(position_length)
            data = bitstring[fragment * data_length:(fragment + 1) * data_length]
            frag = position + data
            encoded = binary2flags(frag, self.flag_count)
            fragments.append(encoded if encoded else 'A')

        return fragments
    
    def decode(self, fragments):
        
        # compute number of bits per fragment
        fragment_length = self.getFragmentLength()
        
        # convert fragments into bitstrings
        fragments = sorted(set(
            flags2binary(fragment, self.flag_count).zfill(fragment_length)
            for fragment in set(fragments)
        ))
        
        # derive number of bits per fragment reserved to encode the position as 
        # one more than the number of leading zeros in the second fragment
        position_length = fragments[1].find('1') + 1
        # position_length = 0
        # while fragments[1][position_length] == '0':
        #     position_length += 1
        # position_length += 1
        
        # concatenate data fragments
        bitstring = ''.join(fragment[position_length:] for fragment in fragments)
        
        # split data fragments in to length bits and data bits
        length, data = bitstring[:self.filesize_length], bitstring[self.filesize_length:]
        
        # remove padding
        return data[:int(length, 2)]