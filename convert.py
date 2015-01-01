#!/usr/bin/env python

from bitarray import bitarray
import sys, getopt, argparse

#define the length in bits of a character
char_length = 8

'''
INPUT:
    character: character as a string
    parity_positions: array of integers representing 0 indexed positions of parity bits for a character
OUTPUT:
    char_hammingArray: bitarray of length char_length plus number of parity bits encoded with (8,12) Hamming Code
'''
def set_hammingArray(character, parity_positions):
    #make a bitarray for the input character
    char_bitarray = bitarray()
    char_bitarray.frombytes(character)
    #make a bit array to hold the channel coded version of the character. its size must be expanded by the number of parity bits
    char_hammingArray = bitarray(char_bitarray.length()+len(parity_positions))

    #transfer all bits from the char_bitarray to the char_hammingArray skipping the parity bit positions
    i = 0
    j = 0
    while (i < char_hammingArray.length()):
        if i in parity_positions:
            i += 1
            continue
        char_hammingArray[i] = char_bitarray[j]
        i += 1
        j += 1

    #set each parity bit for the hamming code: 1 if checksum is odd, 0 if checksum is even
    for bit in parity_positions:
        pos = bit
        bit+=1
        check_sum = 0
        #indicate if skipping or checking while iterating down char_hammingArray
        flip_ind = 0

        #pos starts at the parity bit position each iteration and is incremented by 1 until it exceeds the length of the char_hammingArray
        while (pos < char_hammingArray.length()):
            check_sum += char_hammingArray[pos]
            flip_ind += 1
            if(flip_ind == bit):
                pos += bit
                flip_ind = 0
            pos+=1
        #if the checksum is odd, set the parity bit to 1, otherwise leave at 0
        if(check_sum%2!=0):
            char_hammingArray[(bit-1)] = 1
    return char_hammingArray

'''
Convert 1s and 0s to A,T,C,G with following pattern:
A,T: 0
C,G: 1
bits alternate A,T/C,G to allow for more structurally sound DNA and minor error checking against indels
INPUT: bit array for hamming encoded character
OUTPUT: string of DNA for hamming encoded character
'''
def convert_char_bitarray(char_bitArray):
    #loop through each bit of the hamming code array
    zero_ind = 0
    one_ind = 0
    dna_out = ''
    for bit in char_bitArray:
        if bit:
            if one_ind:
                dna_out += 'C'
                one_ind = 0
            else:
                dna_out += 'G'
                one_ind = 1
        else:
            if zero_ind:
                dna_out += 'A'
                zero_ind = 0
            else:
                dna_out += 'T'
                zero_ind = 1
    return dna_out

def main(argv):
    inputfile = ''
    outputfile = ''
    print "Reading input and output files..."
    #set up and parse command line arguments. Handle errors if not provided properly
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("input")
    parser.add_argument("output", nargs='?', default='out.txt')
    args = parser.parse_args()
    args.input
    inputfile = args.input
    outputfile = args.output

    #open file descriptors of inptu and output files
    input_file_object = open(inputfile, 'r')
    output_file_object = open(outputfile, 'w')
    #convert input file to a string and then split into a list of characters
    input_text = input_file_object.read()
    input_list = list(input_text)
    print "Input File: ",inputfile
    print "Output File: ",outputfile

    print "Converting to Hamming Code..."
    #define an array of positions for parity bits
    parity_positions = []
    parity_bit = 1
    while parity_bit <= char_length:
        parity_positions.append(parity_bit-1)
        parity_bit *= 2

    seq_out = ''
    num_chars = 0
    #for every character in the input list, make a bit array for it encoded with the hamming code
    for character in input_list:
        char_hammingArray = set_hammingArray(character, parity_positions)
        seq_out += convert_char_bitarray(char_hammingArray)
        #every 80 characters, add a new line to adhere to recommended FASTA format
        if num_chars%10==0:
            seq_out+='\n'
        num_chars+=1
    print "Writing output..."
    output_file_object.write('>Hamming Code Output: '+outputfile+', Characters Encoded: '+str(num_chars)+'\n')
    output_file_object.write(seq_out)
    print "Successfully converted"

if __name__ == "__main__":
    main(sys.argv[1:])
