import os, subprocess
import re
# from pathlib import Path
import shutil
from Bio.Align.Applications import MuscleCommandline

class MSA:

    def __init__(self,name,data,folder):
        self.name=name
        self.data=data
        self.folder=folder

    def compute(self):
        output=self.folder+'train_sequences/'
        if self.name == "muscle":
            output=self.compute_muscle()
        elif self.name == "mafft":
            output=self.compute_mafft()
        elif self.name == "tm-coffee":
            output = self.compute_tmcoffee()
        elif self.name == "aqua":
            output=self.compute_aqua()
        elif self.name == "clustalomega":
            output=self.compute_clustalomega()
        elif self.name == "clustalw":
            output=self.compute_clustalw()
        elif self.name == "t_coffee":
            output=self.compute_tcoffee()
        return output

    def compute_aqua(self):
        path_out = self.folder+"Aqua_output/"
        
        
        if os.path.exists(path_out):
            shutil.rmtree(path_out)
            
        os.mkdir(path_out)
        
        blast_fasta = os.walk(self.data, topdown=True)
        list = []
        for path,dir_list,file_list in blast_fasta:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    list.append(file_name)
                    
                else:
                    print ('Invalid File Name :' + file_name)

            for file_name in list:
                input_aqua = path+file_name
                subprocess.call(["singularity", "exec", "-H", "/speed-scratch/shiva/", "/speed-scratch/bioinformatics-group/aqua.simg", "AQUA.tcl",input_aqua, path_out])

        # output=os.walk(path_out, topdown=True)
        # for path,dir_list,file_list in output:
        #     for file_name in file_list:
        #         prog = re.compile('^\d')
        #         result = prog.match(file_name)
        #         if (result):
        #             if file_name[(len(file_name)-4):]=='best':
        #                 new_file=(file_name[:-11])+'.fasta'
        #                 fs = open(file_name, 'r+')
        #                 groups = fs.readlines()
        #                 with open(best+new_file,'a') as newfile:
        #                     for group in groups:
        #                         newfile.writelines(group)
        #                     newfile.flush()

                        # subprocess.call(['mv',path+file_name, best+new_file])
                        # os.rename(path+file_name, best+new_file)
                        # shutil.move(path+file_name, best+new_file)
                        
                    # else:
                        # os.remove(path+file_name)
                        
        # shutil.rmtree(path_out)
        

        best=self.find_best(path_out)
        return(best)

    def find_best(self,pathout):

        best=self.folder+"aqua_result/"
        if os.path.exists(best):
            shutil.rmtree(best)
        os.mkdir(best)
        output=os.walk(pathout, topdown=True)
        for path,dir_list,file_list in output:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    if file_name[(len(file_name)-4):]=='best':
                        new_file=(file_name[:-11])+'.fasta'
                            
                        with open(path+file_name) as f:
                            lines = f.readlines()
                            line = [l for l in lines ]
                            with open(best+new_file, "a") as f1:
                                f1.writelines(line)
        return best


    def compute_clustalw(self):
        path_out = self.folder+"clustalw_result/"
        if os.path.exists(path_out):
            shutil.rmtree(path_out)
        os.mkdir(path_out)
        blast_fasta = os.walk(self.data, topdown=True)
        list = []
        for path, dir_list, file_list in blast_fasta:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    list.append(file_name)
                else:
                    print('Invalid File Name :' + file_name)
            for file_name in list:
                input_clustal = path + file_name
                output_clustal = path_out + file_name[0:(len(file_name) - 6)] + '.fasta'

                subprocess.call(["singularity", "exec", "-H", "/speed-scratch/shiva/", "/speed-scratch/bioinformatics-group/bioinformatics-singularity.simg", "clustalw",'-infile='+ str(input_clustal), '-outfile='+str(output_clustal), '-OUTPUT=FASTA'])

        return (path_out)


    def compute_tcoffee(self):
        path_out = self.folder+"tcoffee_result/"
        if os.path.exists(path_out):
            shutil.rmtree(path_out)

        os.mkdir(path_out)
        blast_fasta = os.walk(self.data, topdown=True)
        list = []
        for path, dir_list, file_list in blast_fasta:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    list.append(file_name)
                else:
                    print('Invalid File Name :' + file_name)
            for file_name in list:
                input_tcoffee = path + file_name
                output_tcoffee = path_out + file_name[0:(len(file_name) - 6)] + '.fasta'

                subprocess.call(["singularity", "exec", "-H", "/speed-scratch/shiva/",
                                 "/speed-scratch/bioinformatics-group/bioinformatics-singularity.simg", "t_coffee",
                                 '-infile', input_tcoffee, '-outfile', output_tcoffee, "-output", "fasta_aln"])

        return (path_out)





    def compute_muscle(self):
        path_out = self.folder+"MUSCLE_result/"
        if os.path.exists(path_out):
            shutil.rmtree(path_out)
        os.mkdir(path_out)
        blast_fasta = os.walk(self.data, topdown=True)
        list = []
        for path,dir_list,file_list in blast_fasta:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    list.append(file_name)
                else:
                    print ('Invalid File Name :' + file_name)
            for file_name in list:
                input_muscle = path + file_name
                output_muscle = path_out + file_name[0:(len(file_name)-6)] + '.fasta'
                with open(output_muscle, 'w') as outputFile:
                    subprocess.call(["singularity", "exec", "-H", "/speed-scratch/shiva/", "/speed-scratch/bioinformatics-group/aqua.simg","muscle", '-in',input_muscle, '-out',output_muscle])

        return(path_out)

    def compute_clustalomega(self):
        path_out = self.folder+"clustalomega_result/"
        if os.path.exists(path_out):
            shutil.rmtree(path_out)
        os.mkdir(path_out)
        blast_fasta = os.walk(self.data, topdown=True)
        list = []
        for path,dir_list,file_list in blast_fasta:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    list.append(file_name)
                else:
                    print ('Invalid File Name :' + file_name)
            for file_name in list:
                input_clustal = path + file_name
                output_clustal = path_out + file_name[0:(len(file_name)-6)] + '.fasta'
                subprocess.call(["singularity", "exec", "-H", "/speed-scratch/shiva/", "/speed-scratch/bioinformatics-group/bioinformatics-singularity.simg","clustalo", '-infile='+ str(input_clustal), '-outfile='+str(output_clustal), '-OUTPUT=FASTA'])

        return(path_out)


    def compute_mafft(self):
        # windows path
        # path_out = str(Path().absolute())+"\MAFFT_result\"
        # mafft_exe = str(Path().absolute())+"\tools\mafft-mac\mafft.bat"
        path_out = self.folder+"MAFFT_result/"
        # mafft_exe = "tools/mafft-mac/mafft.bat"
        if os.path.exists(path_out):
            shutil.rmtree(path_out)
        os.mkdir(path_out)
        c = os.walk(self.data, topdown=True)
        list = []
        for path,dir_list,file_list in c:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    list.append(file_name)
                else:
                    print ('Invalid File Name :' + file_name)
            for file_name in list:
                input_path = path + file_name
                fasta_file_name = file_name[0:(len(file_name)-6)]
                out_put_path = path_out + fasta_file_name + '.fasta'
                with open(out_put_path, 'w') as outputFile:
                    subprocess.call(["singularity", "exec", "-H", "/speed-scratch/shiva/", "/speed-scratch/bioinformatics-group/aqua.simg","mafft",'--auto', input_path], stdout=outputFile)
        return path_out

    def compute_tmcoffee(self):
        # windows path
        # tm_coffee_out = str(Path().absolute())+"\TM-COFFEE_result\"
        tm_coffee_out = self.folder+"TM-COFFEE_result/"
        if os.path.exists(tm_coffee_out):
            shutil.rmtree(tm_coffee_out)
        os.mkdir(tm_coffee_out)
        # windows path
        # make_db_output= str(Path().absolute())+"\database\database"
        make_db_output = self.folder+"database/blast_database"
        d = os.walk(self.data, topdown=True)

        tm_coffee_valid_file_list = []
        for path,dir_list,file_list in d:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    tm_coffee_valid_file_list.append(file_name)
                else:
                    print('Invalid File Name :' + file_name)
            for file_name in tm_coffee_valid_file_list:
                tm_coffee_input_path = path + file_name
                tm_coffee_fasta_file_name = file_name[0:(len(file_name)-6)]
                tm_coffee_out_put_path = tm_coffee_out + tm_coffee_fasta_file_name + '.fasta'
                f= open(tm_coffee_out_put_path,"w+")
                f.close()
                subprocess.call(["singularity", "exec", "-H", "/speed-scratch/shiva/",'/speed-scratch/bioinformatics-group/bioinformatics-singularity.beta.simg',"t_coffee", tm_coffee_input_path, "-outfile", tm_coffee_out_put_path, "-output", "fasta_aln"
                                , "-mode", "psicoffee", "-blast_server=LOCAL", "-protein_db",make_db_output])

        return tm_coffee_out