import os,subprocess
import shutil
from Bio import SeqIO

class DATABASE:

    def __init__(self, trainset, testset, dataset,folder):
        self.trainset=trainset
        self.testset = testset
        self.dataset = dataset
        self.folder = folder

    def clean(self):
        data_set={}
        train_test = {}

        for seq_record in SeqIO.parse(self.dataset, 'fasta'):
            data_set[seq_record.id]=seq_record.seq

        for seq_record in SeqIO.parse(self.trainset, 'fasta'):
            train_test[seq_record.id]=seq_record.seq
        for seq_record in SeqIO.parse(self.testset, 'fasta'):
            train_test[seq_record.id]=seq_record.seq

        for i in data_set:
            if i in train_test:
                del data_set[i]

        clean_dataset = self.folder + "database/blast_database"
        if os.path.exists(self.folder + "database/"):
            shutil.rmtree(self.folder + 'database/')
        os.makedirs(self.folder + "database/")
        with open(clean_dataset, "w") as result:
            for i in data_set.keys():
                result.write('>'+i+'\n')
                result.write(str(data_set[i]))
        return clean_dataset

    def makedb(self, input):
        make_db_output= self.folder+"database/blast_database"
        subprocess.call(["singularity", "exec", "-H", "/speed-scratch/shiva/",'/speed-scratch/bioinformatics-group/bioinformatics-singularity.simg',"makeblastdb","-in",input,"-title",make_db_output,"-dbtype","prot"])

    def compute(self):
        output=self.clean()
        self.makedb(output)
        os.remove(output)

