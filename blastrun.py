import subprocess
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from Bio import SearchIO
import re
from os.path import isfile, join
import os
import shutil
import sys
# from pathlib import Path

class BlastRun:

    def __init__(self,folder):
        self.folder=folder

        self.blast_run = self.folder+"blastxml/"
        self.blast_fasta_folder_uncheck = self.folder+"blastfasta_uncheck/"
        self.blast_fasta_folder=self.folder+'blastfasta/'
        self.queries = self.folder+"train_sequences/"
        self.make_db_output = self.folder+"database/blast_database"

    def run(self):
        

        blast_out_path = self.blast_run
        a = os.walk(self.queries)
        blast_valid_file_list_first = []
        for path,dir_list,file_list in a:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    blast_valid_file_list_first.append(file_name)
                else:
                    print ('Invalid File Name :' + file_name)
            for file_name in blast_valid_file_list_first:
                blast_query_path = path + file_name
                blast_xml = blast_out_path + file_name[0:(len(file_name)-6)] + '.xml'
                # blastp_cline = NcbiblastpCommandline('blastp', query=blast_query_path, db=make_db_output, evalue=0.001,
                #                                      outfmt=5, out=blast_xml, max_target_seqs=10)
                subprocess.call(["singularity", "exec", "-H", "/speed-scratch/shiva/", '/speed-scratch/bioinformatics-group/bioinformatics-singularity.simg','blastp','-query',blast_query_path, '-db',self.make_db_output, '-evalue','0.001',
                                                     '-outfmt','5', '-out',blast_xml, '-max_target_seqs','11'],stdout=subprocess.PIPE)
                # subprocess.call(str(blastp_cline), stdout=subprocess.PIPE, shell=True)

    def xml(self):
        xml_files = self.blast_run
        b = os.walk(xml_files)
        blast_valid_file_list = []
        for path, dir_list, file_list in b:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if result:
                    blast_valid_file_list.append(file_name)
                else:
                    print('Invalid File Name :' + file_name)
            for file_name in blast_valid_file_list:
                blast_xml_file = path + file_name
                blast_fasta_file = self.blast_fasta_folder_uncheck + file_name[0:(len(file_name) - 4)] + '.fasta'
                blast_qresult = SearchIO.read(blast_xml_file, 'blast-xml')
                records = []
                for hit in blast_qresult:
                    records.append(hit[0].hit)
                SeqIO.write(records, blast_fasta_file, "fasta")
                del records[:]

    def check(self):

        blast_out_path = self.blast_fasta_folder_uncheck
        a = os.walk(blast_out_path)
        
        for path,dir_list,file_list in a:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                count=-1
                for seq_record in (SeqIO.parse(blast_out_path+file_name, 'fasta')):
                    count+=1
                if count!=0:
                            
                    with open(self.blast_fasta_folder_uncheck+file_name) as f:
                        lines = f.readlines()
                        line = [l for l in lines ]
                        with open(self.blast_fasta_folder+file_name, "a") as f1:
                            f1.writelines(line)
                    

    def compute(self):
        if os.path.exists(self.blast_run):
            shutil.rmtree(self.blast_run)

        if os.path.exists(self.blast_fasta_folder_uncheck):
            shutil.rmtree(self.blast_fasta_folder_uncheck)
        
        if os.path.exists(self.blast_fasta_folder):
            shutil.rmtree(self.blast_fasta_folder)

        os.makedirs(self.blast_run)
        os.makedirs(self.blast_fasta_folder_uncheck)
        os.makedirs(self.blast_fasta_folder)

        self.run()
        self.xml()
        self.check()
        return(self.blast_fasta_folder)




