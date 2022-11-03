# On my honor, I have neither given nor received unauthorized aid on this assignment

# Prepared by : Arman Singhal
# UFID: 6602-1080
# Email: armansinghal@ufl.edu

# Import Statements
from bitarray import bitarray
from bitarray.util import ba2int
import numpy as np
import sys

# File Reading
file = open(sys.argv[1], 'r')
content = file.readlines()

regsiters = [0] * 32
data = []

disassembly = {}

# Classifying Instructions and Data Values
def set_type(line):
    bin_inst = disassembly[line]["binary_instruct"]
    if bin_inst[0:6] == "000110":
        for i in range(line):
            disassembly[i+1]["line_type"] = "Instruction"
        
        for i in range(line, len(disassembly)):
            disassembly[i+1]["line_type"] = "Data"

# Category 1 type Instruction Function
def category1_disassembler(line):
    opcode_bits = disassembly[line]["binary_instruct"][3:6]
    rs = disassembly[line]["binary_instruct"][6:11]
    rt = disassembly[line]["binary_instruct"][11:16]
    offset = str(disassembly[line]["binary_instruct"][16:]) + '00'
    rs_ins = "R" + str(int(rs, 2))
    rt_ins = "R" + str(int(rt, 2))
    dest = "#" + str(int(offset, 2))
    instruction = ""
    
    if opcode_bits == '000':
        instruction = instruction + "J "
        target_address = str(disassembly[line]["binary_instruct"][6:]) + '00'
        dest = "#" + str(int(target_address, 2))
        instruction = instruction + dest
    elif opcode_bits == '001':
        instruction = instruction + "BEQ "
        instruction = instruction + rs_ins + ", " + rt_ins + ", " + dest
    elif opcode_bits == '010':
        instruction = instruction + "BNE "
        instruction = instruction + rs_ins + ", " + rt_ins + ", " + dest
    elif opcode_bits == '011':
        instruction = instruction + "BGTZ "
        instruction = instruction + rs_ins + ", " + dest
    elif opcode_bits == '100':
        instruction = instruction + "SW "
        offset = str(disassembly[line]["binary_instruct"][16:])
        dest = str(int(offset, 2))
        instruction = instruction + rt_ins + ", " + dest + "(" + rs_ins + ")"
    elif opcode_bits == '101':
        instruction = instruction + "LW "
        offset = str(disassembly[line]["binary_instruct"][16:])
        dest = str(int(offset, 2))
        instruction = instruction + rt_ins + ", " + dest + "(" + rs_ins + ")"
    elif opcode_bits == '110':
        instruction = instruction + "BREAK "

    return instruction

# Category 2 type Instruction Function
def category2_disassembler(line):
    opcode_bits = disassembly[line]["binary_instruct"][3:6]
    dest_bits = disassembly[line]["binary_instruct"][6:11]
    src1_bits = disassembly[line]["binary_instruct"][11:16]
    src2_bits = disassembly[line]["binary_instruct"][16:21]
    b = bitarray(src2_bits)
    data_val = ba2int(b, signed=True)
    instruction = ""
    
    if opcode_bits == '000':
        instruction = instruction + "ADD "
    elif opcode_bits == '001':
        instruction = instruction + "SUB "
    elif opcode_bits == '010':
        instruction = instruction + "AND "
    elif opcode_bits == '011':
        instruction = instruction + "OR "
    elif opcode_bits == '100':
        instruction = instruction + "SRL "
    elif opcode_bits == '101':
        instruction = instruction + "SRA "
    elif opcode_bits == '110':
        instruction = instruction + "MUL "
    
    dest = "R" + str(int(dest_bits, 2)) + ", "
    instruction = instruction + dest
    
    src1 = "R" + str(int(src1_bits, 2)) + ", "
    instruction = instruction + src1
    
    if opcode_bits == '100' or opcode_bits == '101':
        src2 = "#" + str(data_val)
    else:
        src2 = "R" + str(data_val)
    instruction = instruction + src2

    return instruction

# Category 3 type Instruction Function
def category3_disassembler(line):
    opcode_bits = disassembly[line]["binary_instruct"][3:6]
    instruction = ""
    dest_bits = disassembly[line]["binary_instruct"][6:11]
    src_bits = disassembly[line]["binary_instruct"][11:16]
    immVal_bits = disassembly[line]["binary_instruct"][16:33]
    b = bitarray(immVal_bits)
    data_val = ba2int(b, signed=True)
    
    if opcode_bits == '000':
        instruction = instruction + "ADDI "
    elif opcode_bits == '001':
        instruction = instruction + "ANDI "
    elif opcode_bits == '010':
        instruction = instruction + "ORI "
    
    dest = "R" + str(int(dest_bits, 2)) + ", "
    instruction = instruction + dest
    
    src = "R" + str(int(src_bits, 2)) + ", "
    instruction = instruction + src
    
    immVal = "#" + str(data_val)
    instruction = instruction + immVal
        
    return instruction

# Disassebler Routine
disassembly = {}
for i in range(len(content)):
    disassembly[i+1] = {"binary_instruct":content[i][0:32], "mem_address":((i * 4) + 260)}

for line in disassembly:
    set_type(line)

data_dict = {}
for line in disassembly:
    if disassembly[line]["line_type"] == "Data":
        data_dict[line] = disassembly[line]

disassembly = {k: v for k, v in disassembly.items() if k not in data_dict}

with open("disassembly.txt", "w") as f:
    for line in disassembly:
        cat_msb = disassembly[line]["binary_instruct"][0:3]
        if cat_msb == '000':
            disassembly[line]["instruction"] = category1_disassembler(line)
        elif cat_msb == '001':
            disassembly[line]["instruction"] = category2_disassembler(line)
        elif cat_msb == "010":
            disassembly[line]["instruction"] = category3_disassembler(line)
        write_instr = str(disassembly[line]["binary_instruct"]) + "\t" + str(disassembly[line]["mem_address"]) + "\t" + disassembly[line]["instruction"]
        f.write(write_instr)
        f.write("\n")

    for line in data_dict:
        b = bitarray(data_dict[line]["binary_instruct"])
        data_val = ba2int(b, signed=True)
        data_dict[line]["value"] = data_val
        write_line = str(data_dict[line]["binary_instruct"]) + "\t" + str(data_dict[line]["mem_address"]) + "\t" + str(data_val)
        f.write(write_line)
        f.write("\n")

for line in data_dict:
    data.append(data_dict[line]["value"])

# Simulator Functions

# Tested and Approved
def add_inst(dest, src1, src2):
    dest = int(dest[1:-1])
    src1 = int(src1[1:-1])
    src2 = int(src2.split("R")[1])
    regsiters[dest] = regsiters[src1] + regsiters[src2]

# Tested and Approved
def sub_inst(dest, src1, src2):
    dest = int(dest[1:-1])
    src1 = int(src1[1:-1])
    src2 = int(src2.split("R")[1])
    regsiters[dest] = regsiters[src1] - regsiters[src2]

def and_inst(dest, src1, src2):
    dest = int(dest[1:-1])
    src1 = int(src1[1:-1])
    src2 = int(src2.split("R")[1])
    regsiters[dest] = regsiters[src1] & regsiters[src2]

def or_inst(dest, src1, src2):
    dest = int(dest[1:-1])
    src1 = int(src1[1:-1])
    src2 = int(src2.split("R")[1])
    regsiters[dest] = regsiters[src1] | regsiters[src2]

def srl_inst(dest, src1, src2):
    dest = int(dest[1:-1])
    src1 = int(src1[1:-1])
    src2 = int(src2[1:])
    regsiters[dest] = np.right_shift(regsiters[src1], src2)

def sra_inst(dest, src1, src2):
    dest = int(dest[1:-1])
    src1 = int(src1[1:-1])
    src2 = int(src2[1:])
    regsiters[dest] = regsiters[src1] >> src2

# Tested and Approved
def mul_inst(dest, src1, src2):
    dest = int(dest[1:-1])
    src1 = int(src1[1:-1])
    src2 = int(src2[1:])
    regsiters[dest] = regsiters[src1] * regsiters[src2]

# Tested and Approved
def j_inst(inst_address):
    mem = int(inst_address[1:])
    return int(((mem-260)/4)+1) - 1

# Tested and Approved
def beq_inst(rs, rt, offset, mem_):
    branch_mem = mem_
    rs = int(rs[1:-1])
    rt = int(rt[1:-1])
    offset = int(offset[1:])
    if regsiters[rs] == regsiters[rt]:
        branch_mem = mem_ + 4 + offset
        return int((branch_mem-260)/4)
    else:
        return int((branch_mem-260+4)/4)

def bne_inst(rs, rt, offset, mem_):
    branch_mem = mem_
    rs = int(rs[1:-1])
    rt = int(rt[1:-1])
    offset = int(offset[1:])
    if regsiters[rs] != regsiters[rt]:
        branch_mem = mem_ + 4 + offset
        return int((branch_mem-260)/4)
    else:
        return int((branch_mem-260+4)/4)

# Tested and Approved
def addi_inst(dest, src1, src2):
    dest = int(dest[1:-1])
    src1 = int(src1[1:-1])
    src2 = int(src2[1:])
    regsiters[dest] = regsiters[src1] + src2

# Tested and Approved
def sw_inst(src, offset):
    src = int(src[1:-1])
    base = int(offset.split("R")[1][:-1])
    mem_ = int(offset.split("R")[0][:-1]) + regsiters[base]
    first_data_mem = list(data_dict.values())[0]["mem_address"]
    data[int((mem_ - first_data_mem)/4)] = regsiters[src]

# Tested and Approved
def lw_inst(dest, offset):
    dest = int(dest[1:-1])
    base = int(offset.split("R")[1][:-1])
    mem_ = int(offset.split("R")[0][:-1]) + regsiters[base]
    first_data_mem = list(data_dict.values())[0]["mem_address"]
    load_data = data[int((mem_ - first_data_mem)/4)]
    regsiters[dest] = load_data

# Tested and Approved
def bgtz_inst(src, offset, mem_):
    branch_mem = mem_
    src = int(src[1:-1])
    offset = int(offset[1:])
    if regsiters[src] > 0:
        branch_mem = mem_ + 4 + offset
        return int((branch_mem-260)/4)
    else:
        return int((branch_mem-260+4)/4)

def andi_inst(rt, rs, imm):
    rt = int(rt[1:-1])
    rs = int(rs[1:-1])
    imm = int(imm[1:]).astype('Int32')
    regsiters[rt] = regsiters[rs] & imm


def ori_inst(rt, rs, imm):
    rt = int(rt[1:-1])
    rs = int(rs[1:-1])
    imm = int(imm[1:]).astype('Int32')
    regsiters[rt] = regsiters[rs] | imm

# Simulator Main
cycle = 1
line = 1
with open("simulation.txt", "w") as f:
    while line < len(disassembly)+1:
        mem_add = disassembly[line]["mem_address"]
        instruct = disassembly[line]["instruction"]

        inst_parse = instruct.split(" ")
        keyword = inst_parse[0]

        if keyword == "J":
            line = j_inst(inst_parse[1])
        elif keyword == "BEQ":
            line = beq_inst(inst_parse[1], inst_parse[2], inst_parse[3], mem_add)
        elif keyword == "BNE":
            line = bne_inst(inst_parse[1], inst_parse[2], inst_parse[3], mem_add)
        elif keyword == "LW":
            lw_inst(inst_parse[1], inst_parse[2])
        elif keyword == "SW":
            sw_inst(inst_parse[1], inst_parse[2])
        elif keyword == "BGTZ":
            line = bgtz_inst(inst_parse[1], inst_parse[2], mem_add)
        elif keyword == "ADD":
            add_inst(inst_parse[1], inst_parse[2], inst_parse[3])
        elif keyword == "SUB":
            sub_inst(inst_parse[1], inst_parse[2], inst_parse[3])
        elif keyword == "AND":
            and_inst(inst_parse[1], inst_parse[2], inst_parse[3])
        elif keyword == "OR":
            or_inst(inst_parse[1], inst_parse[2], inst_parse[3])
        elif keyword == "SRL":
            srl_inst(inst_parse[1], inst_parse[2], inst_parse[3])
        elif keyword == "SRA":
            sra_inst(inst_parse[1], inst_parse[2], inst_parse[3])
        elif keyword == "MUL":
            mul_inst(inst_parse[1], inst_parse[2], inst_parse[3])
        elif keyword == "ADDI":
            addi_inst(inst_parse[1], inst_parse[2], inst_parse[3])
        elif keyword == "ANDI":
            andi_inst(inst_parse[1], inst_parse[2], inst_parse[3])
        elif keyword == "ORI":
            ori_inst(inst_parse[1], inst_parse[2], inst_parse[3])


        f.write("--------------------\n")
        f.write("Cycle " + str(cycle) + ":\t" + str(mem_add) + "\t" + instruct + "\n\n")
        f.write("Registers\n")
        reg1 = "R00:\t" + str(regsiters[0]) + "\t" + str(regsiters[1]) + "\t" + str(regsiters[2]) + "\t" + str(regsiters[3]) + "\t" + str(regsiters[4]) + "\t" + str(regsiters[5]) + "\t" + str(regsiters[6]) + "\t" + str(regsiters[7]) + "\n"
        f.write(reg1)
        
        reg2 = "R08:\t" + str(regsiters[8]) + "\t" + str(regsiters[9]) + "\t" + str(regsiters[10]) + "\t" + str(regsiters[11]) + "\t" + str(regsiters[12]) + "\t" + str(regsiters[13]) + "\t" + str(regsiters[14]) + "\t" + str(regsiters[15]) + "\n"
        f.write(reg2)
        
        reg3 = "R16:\t" + str(regsiters[16]) + "\t" + str(regsiters[17]) + "\t" + str(regsiters[18]) + "\t" + str(regsiters[19]) + "\t" + str(regsiters[20]) + "\t" + str(regsiters[21]) + "\t" + str(regsiters[22]) + "\t" + str(regsiters[23]) + "\n"
        f.write(reg3)
        
        reg4 = "R24:\t" + str(regsiters[24]) + "\t" + str(regsiters[25]) + "\t" + str(regsiters[26]) + "\t" + str(regsiters[27]) + "\t" + str(regsiters[28]) + "\t" + str(regsiters[29]) + "\t" + str(regsiters[30]) + "\t" + str(regsiters[31]) + "\n"
        f.write(reg4)
        
        f.write("\nData\n")
        
        first_data_mem = list(data_dict.values())[0]["mem_address"]
        first_mem = str(first_data_mem) + ": "
        for i in range(8):
            first_mem = first_mem + "\t" + str(data[i])
        f.write(first_mem)
        f.write("\n")
        
        no_of_data_values = len(data_dict)
        last_data_index = 8
        sec_data = str(first_data_mem + 32) + ": "
        for i in range(int(no_of_data_values/8) - 1):
            sec_data = sec_data + "\t" + str(data[last_data_index]) + "\t" + str(data[last_data_index+1]) + "\t" + str(data[last_data_index+2]) + "\t" + str(data[last_data_index+3]) + "\t" + str(data[last_data_index+4]) + "\t" + str(data[last_data_index+5]) + "\t" + str(data[last_data_index+6]) + "\t" + str(data[last_data_index+7])
            f.write(sec_data)
            last_data_index += 8
        f.write("\n\n")
        cycle += 1
        line = line + 1