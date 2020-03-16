import os

if __name__ == "__main__":
    for file in os.listdir():
        if file.startswith("test"):
            os.system("python3 {}".format(file))
    print("All tests passed :)")