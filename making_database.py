import os,subprocess
import shutil

class DATABASE:

    def __init__(self, testset, dataset,folder):
        self.testset = testset
        self.dataset = dataset
        self.folder = folder

    def removedup(self):
        small_content = []
        larget_content = {}
        with open(self.testset, "r") as smallFile:
            oneline = smallFile.readline()
            temp = []
            while oneline:
                if not oneline.strip("\n").strip():
                    oneline = smallFile.readline()
                    continue
                if oneline.startswith(">"):
                    if temp:
                        small_content.append("".join(temp))
                        temp = []
                else:
                    temp.append(oneline.strip('\n'))
                oneline = smallFile.readline()
        small_content.append("".join(temp))
        with open(self.dataset, "r") as largeFile:
            oneline = largeFile.readline()
            temp = []
            entire = []
            while oneline:
                if not oneline.strip("\n").strip():
                    oneline = largeFile.readline()
                    continue
                if oneline.startswith(">"):
                    if entire:
                        larget_content.update({"".join(temp): "".join(entire)})
                        temp = []
                        entire = []
                else:
                    temp.append(oneline.strip("\n"))
                entire.append(oneline)
                oneline = largeFile.readline()
        larget_content.update({"".join(temp): "".join(entire)})

        # print("Before removing, the number of larger dataset is "+ str(len(larget_content)))

        for each in small_content:
            if each in larget_content:
                del larget_content[each]

        # print("After removing duplicate records, the number of larger dataset is " + str(len(larget_content)))
        clean_dataset = self.folder+"database/blast_database"
        if os.path.exists(self.folder+"database/"):
            shutil.rmtree(self.folder+'database/')
        os.makedirs(self.folder+"database/")
        with open(clean_dataset, "w") as result:
            result.writelines(larget_content.values())

        return clean_dataset

    def makedb(self, input):
        make_db_output= self.folder+"database/blast_database"
        subprocess.call(["singularity", "exec", "-H", "/speed-scratch/shiva/",'/speed-scratch/bioinformatics-group/bioinformatics-singularity.simg',"makeblastdb","-in",input,"-title",make_db_output,"-dbtype","prot"])

    def compute(self):
        output = self.removedup()
        self.makedb(output)
        os.remove(output)

