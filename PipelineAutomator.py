import subprocess
import time

time1 = time.time()

files = [
    "/home/arvin/barseq/controls/Control_1.mat.h5ad",
    "/home/arvin/barseq/controls/Control_2.mat.h5ad",
    "/home/arvin/barseq/controls/Control_3.mat.h5ad",
    "/home/arvin/barseq/controls/Control_4.mat.h5ad",
    "/home/arvin/barseq/enucleated/Enucleated_1.mat.h5ad",
    "/home/arvin/barseq/enucleated/Enucleated_2.mat.h5ad",
    "/home/arvin/barseq/enucleated/Enucleated_3.mat.h5ad",
    "/home/arvin/barseq/enucleated/Enucleated_4.mat.h5ad"]

for i in range(len(files) - 1):
    file1 = files[i]
    file2 = files[i + 1]
    print(f"Matching {file1} with {file2}")
    subprocess.run(["python", "MatchingByRegion.py", file1, file2])

subprocess.run(["python", "MatchAggregator.py"])

print(f"All Files Done. Took {time.time() - time1}")
