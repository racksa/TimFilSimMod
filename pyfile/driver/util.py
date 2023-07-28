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