import os
import glob

# Missing tests for the C wrapper of SLICOT routines

# Find all *.cpp files in the tests directory
# loop through all cpp files and split the file name at _ and take the first part and append it to a set
# Find all *.c file src_c_wrapper directory and createa set of the file names without the .c extension
# Compare the two sets and print the missing files

def get_tested_routines(tests_dir):
    tested = set()
    for filepath in glob.glob(os.path.join(tests_dir, "*.cpp")):
        filename = os.path.basename(filepath)
        routine = filename.split('_')[0]
        tested.add(routine)
    return tested

def get_c_wrappers(src_dir):
    wrappers = set()
    for filepath in glob.glob(os.path.join(src_dir, "*.c")):
        filename = os.path.basename(filepath)
        routine = os.path.splitext(filename)[0]
        wrappers.add(routine)
    return wrappers

def main():
    tests_dir = os.path.join(os.path.dirname(__file__), "..", "tests")
    src_dir = os.path.join(os.path.dirname(__file__), "..", "src_c_wrapper")

    tested = get_tested_routines(tests_dir)
    wrappers = get_c_wrappers(src_dir)

    missing = wrappers - tested

    print("Missing tests for the following C wrappers:")
    for routine in sorted(missing):
        print(routine)

if __name__ == "__main__":
    main()
