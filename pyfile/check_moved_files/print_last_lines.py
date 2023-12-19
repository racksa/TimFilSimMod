# def print_first_and_last_lines(file_name):
#     line = 0
#     file = open(file_name, 'r')
#     length = sum(1 for line in open(file_name))
#     print(length)
            
#     for i in range(length):
#         l = file.readline().split()

#         if i <= 0:
#             for j in range(10):
#                 print(l[j])
#         if i == length -1:
#             for j in range(10):
#                 print(l[j])

# # Usage
# file_name = "/home/clustor2/ma/h/hs2216/20231009/ciliate_4096fil_10000blob_20.00R_2.00torsion_blob_forces.dat"
# file_name = "ciliate_4096fil_10000blob_20.00R_2.00torsion_blob_forces.dat"
# print_first_and_last_lines(file_name)

def print_last_10_lines(file_name):
    try:
        with open(file_name, 'r') as file:
            file_size = os.path.getsize(file_name)
            last_10_lines = []
            line_count = 0

            # Start from the end of the file and read lines in reverse order
            file.seek(0, os.SEEK_END)
            pointer = file.tell()
            while pointer >= 0 and line_count < 10:
                file.seek(pointer)
                buffer = file.read(1)
                if buffer == "\n":
                    last_10_lines.append(file.readline())
                    line_count += 1
                pointer -= 1

            # Reverse the list of lines to print them in correct order
            last_10_lines = last_10_lines[::-1]

            print("Last 10 Lines:")
            print("".join(last_10_lines))
    except FileNotFoundError:
        print(f"File '{file_name}' not found.")

# Usage
import os
# file_name = "ciliate_4096fil_10000blob_20.00R_2.00torsion_blob_forces.dat"  # Replace with the name of your file
file_name = "/home/clustor2/ma/h/hs2216/20231009/ciliate_4096fil_10000blob_20.00R_2.00torsion_blob_forces.dat"
print_last_10_lines(file_name)
