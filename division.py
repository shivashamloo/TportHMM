from Bio import SeqIO
import os
import shutil
# from pathlib import Path


class DIVISION:

    def __init__(self,data, type,folder):
        self.data = data
        self.type = type
        self.folder=folder

    def divide(self):

        file_type = "fasta"
        filetoseq={}
        
        if(self.type == "train"):
            # windows path
            # original_output = str(Path().absolute())+"\train_sequences\"
            original_output = self.folder+"train_sequences/"
            mapping=self.folder+"trainmap.txt"
        elif (self.type == "test"):
            # windows path
            # original_output = str(Path().absolute())+"\test_sequences\"
            original_output = self.folder+"test_sequences/"
            mapping=self.folder+"testmap.txt"
        if os.path.exists(original_output):
            shutil.rmtree(original_output)
        os.makedirs(original_output)



        for f, seq_record in enumerate(SeqIO.parse(self.data, file_type)):
            o_output = original_output + str(f) + ".fasta"
            SeqIO.write(seq_record, o_output, "fasta")
            filetoseq[f+1]=seq_record.id

        with open(mapping, "w") as f:
            for key, value in filetoseq.items():
                f.write('%s:%s\n' % (key, value))
        

