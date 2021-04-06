import os
import sys
import subprocess
import platform


class TestCpp:
    def filter_executables(self):
        plt = platform.system()
        dir_with_files = os.getcwd()
        if plt == "Linux":
            list_test = [f for f in sorted(os.listdir(dir_with_files)) if (str(f))[-1:] == 'x']
            self.execute_executables(list_test)
        elif plt == "Windows":
            # print("win")
            list_test = [f for f in sorted(os.listdir(dir_with_files)) if (str(f))[-3:] == 'exe']
            self.execute_executables(list_test)

    @staticmethod
    def execute_executables(list_test):
        f = open("testcpp_output.txt", "w")
        for i in range(len(list_test)):
            if list_test[i] == "test_libFLAME.exe":
                continue

            if platform.system() == "Linux":
                command = "./" + list_test[i]
            elif platform.system() == "Windows":
                command = list_test[i]

            print("\n" + "=" * 40 + " " + list_test[i] + " execution " + "=" * 40 + "\n")
            f.write("\n" + "=" * 40 + " " + list_test[i] + " execution " + "=" * 40 + "\n")
            process = subprocess.Popen(command, bufsize=1, universal_newlines=True,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT)
            # file_name = list_test[i] + "_output.txt"

            for line in iter(process.stdout.readline, ''):
                print(line[:-1])
                f.write(line)
                sys.stdout.flush()
            process.wait()
            errcode = process.returncode
        f.close()


if __name__ == '__main__':
    test_object = TestCpp()
    test_object.filter_executables()
