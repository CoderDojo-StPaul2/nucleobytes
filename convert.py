#!/usr/bin/env python

from bitarray import bitarray
import sys, getopt, argparse



def main(argv):
    inputfile = ''
    outputfile = ''
    #define the length in bits of a character
    char_length = 8

    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("input")
    parser.add_argument("output", nargs='?', default='out.txt')
    args = parser.parse_args()
    args.input
    inputfile = args.input
    outputfile = args.output
    print "Input File: ",inputfile
    print "Output File: ",outputfile

    input_file_object = open(inputfile, 'r')
    output_file_object = open(outputfile, 'w')

    input_text = input_file_object.read()
    input_list = list(input_text)
    print input_list
    #define an array of positions for parity bits
    parity_positions = []
    parity_bit = 1
    while parity_bit <= char_length:
        parity_positions.append(parity_bit-1)
        parity_bit *= 2

    for character in input_list:
        char_bitarray = bitarray()
        char_bitarray.frombytes(character)
        #make a bit array to hold the channel coded version. its size must be expanded by the number of parity bits
        char_hammingArray = bitarray(char_bitarray.length()+len(parity_positions))
        i = 1

        while (i <= char_bitarray.length()):
            char_hammingArray[(i-1)] = 1
            i = i*2

        i = 0
        j = 0
        while (i < char_hammingArray.length()):
            if char_hammingArray[i]:
                char_hammingArray[i] = 0
                i += 1
                continue
            char_hammingArray[i] = char_bitarray[j]
            i += 1
            j += 1

        for bit in parity_positions:
            pos = bit
            bit+=1
            bit_sum = 0
            flip_ind = 0

            while (pos < char_hammingArray.length()):
                bit_sum += char_hammingArray[pos]
                flip_ind += 1
                if(flip_ind == bit):
                    pos += bit
                    flip_ind = 0
                pos+=1
            if(bit_sum%2!=0):
                char_hammingArray[(bit-1)] = 1


        print char_hammingArray
        print char_bitarray

if __name__ == "__main__":
    main(sys.argv[1:])
