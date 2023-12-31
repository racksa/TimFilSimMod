import os

def list_files_in_directory(directory_path):
    file_list = []
    for file in os.listdir(directory_path):
        file_path = os.path.join(directory_path, file)
        if os.path.isfile(file_path):
            file_list.append(file_path)
    return file_list


def delete_files_in_directory(directory_path):
    try:
        file_list = list_files_in_directory(directory_path)

        if not file_list:
            print(f"No files found in '{directory_path}'. Nothing to delete.")
            return

        print("Files to be deleted:")
        for file_path in file_list:
            print(file_path)

        user_input = input("Do you want to delete these files? (y/n): ")
        if user_input.lower() == 'y':
            for file_path in file_list:
                os.remove(file_path)
            print("All files have been deleted.")
        else:
            print("Deletion canceled. No files were deleted.")
    except Exception as e:
        print(f"Error occurred while deleting files: {e}")

def view_files_in_directory(directory_path):
    try:
        file_list = list_files_in_directory(directory_path)

        if not file_list:
            print(f"No files found in '{directory_path}'.")
            return

        print("Files here:")
        for file_path in file_list:
            print(file_path)

    except Exception as e:
        print(f"Error occurred while viewing files: {e}")

def copy_last_line(input_filename, output_filename):
    try:
        # Open the input file in read mode
        with open(input_filename, 'r') as input_file:
            # Read all lines from the file
            lines = input_file.readlines()

            # Check if the file is not empty
            if lines:
                # Extract the last line
                last_line = lines[-1]

                # Open the output file in write mode
                with open(output_filename, 'w') as output_file:
                    # Write the last line to the output file
                    output_file.write(last_line)

                print(f"Last line copied from '{input_filename}' to '{output_filename}'.")
            else:
                print(f"The file '{input_filename}' is empty.")
    except FileNotFoundError:
        print(f"Error: The file '{input_filename}' does not exist.")


def even_list_index(n, m):
    sublist_length = n // m  # Floor division to get the length of each sublist
    remainder = n % m  # Get the remainder to distribute extra elements

    result = [0]
    start = 0

    for i in range(m):
        end = start + sublist_length + (1 if i < remainder else 0)
        result.append(end)
        start = end

    return result