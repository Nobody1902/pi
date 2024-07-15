file01 = input("Enter correct pi file: ")
file02 = input("Enter test pi file: ")

pi01 = ""
pi02 = ""

with open(file01, "r") as f:
    pi01 = f.read()

with open(file02, "r") as f:
    pi02 = f.read()

size = min(len(pi01), len(pi02))

for i in range(size):
    if pi01[i] != pi02[i]:
        print(f"Error found at {i} char.")
        print(f"Found '{pi02[i]}' instead of '{pi01[i]}'.")
        exit()

print(f"No error found in {size} digits.")
