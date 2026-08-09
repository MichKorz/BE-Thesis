import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: python replace_commas.py <filename>")
        sys.exit(1)

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            content = file.read()

        modified_content = content.replace(',', '.')

        with open(file_path, 'w', encoding='utf-8') as file:
            file.write(modified_content)

        print(f"Successfully replaced all commas in '{file_path}'.")

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()