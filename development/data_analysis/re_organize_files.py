import shutil
import pathlib
import os


def get_depth(path: str | pathlib.Path, src_dir: str | pathlib.Path) -> int:
    path = str(path)
    src_dir = str(src_dir)
    if path == src_dir:
        return 0
    return len(os.path.relpath(path, src_dir).split(os.sep)) + 1


def copy_and_rename_files(src_dir, dst_dir, files: list[str], max_depth: int = 0):
    # Create destination directory if it doesn't exist
    os.makedirs(dst_dir, exist_ok=True)

    # Walk through the source directory
    for root, dirs, _ in os.walk(src_dir):
        if get_depth(root, src_dir) > max_depth:
            break

        for dir_name in dirs:
            subfolder_path = os.path.join(root, dir_name)
            paths = [os.path.join(subfolder_path, file) for file in files]

            if all(os.path.exists(path) for path in paths):
                dir_name_clean = dir_name.replace(".D", "")
                names = [f"{dir_name_clean}_{file}" for file in files]
                new_paths = [os.path.join(dst_dir, name) for name in names]

                # Copy the files to the new destination with the new names
                for old_path, new_path in zip(paths, new_paths):
                    shutil.copy2(old_path, new_path)
                print(f"Copied {dir_name} to {dst_dir}")
            else:
                print(f"Skipped {dir_name} in {subfolder_path} as not files found.")


def main():
    # Usage
    src_directory = pathlib.Path(r"D:\11_52")
    dst_directory = pathlib.Path(r"C:\Users\nicep\Desktop\research_wis\data\11\11_52\GCMS")
    print(dst_directory)
    files = ['data.ms', 'FID1A.ch', 'pre_post.ini']
    copy_and_rename_files(src_directory, dst_directory, files)


if __name__ == '__main__':
    main()
