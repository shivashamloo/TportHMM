import os, subprocess,sys
import re
# from pathlib import Path
import shutil
from os.path import isfile, join
from Bio import SeqIO


class SDS:

    def __init__(self,name,data,folder):
        self.name = name
        self.data = data
        self.folder=folder

    def compute(self):
        output=self.data
        if self.name == "tcs":
            output = self.compute_tcs()
        elif self.name == "xdet":
            output = self.compute_xdet()
        elif self.name == "speerserver":
            output = self.compute_speer()
        elif self.name == "groupsim":
            output = self.compute_groupsim()
        return output

    def compute_tcs(self):
        
        tm_coffee_out = self.folder+"TCS_result/"
        if os.path.exists(tm_coffee_out):
            shutil.rmtree(tm_coffee_out)

        os.mkdir(tm_coffee_out)
        
        d = os.walk(self.data, topdown=True)
        tm_coffee_valid_file_list = []
        for path, dir_list, file_list in d:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if result:
                    tm_coffee_valid_file_list.append(file_name)
                else:
                    print('Invalid File Name :' + file_name)
            for file_name in tm_coffee_valid_file_list:
                tm_coffee_input_path = path + file_name
                
                tm_coffee_out_put_path = tm_coffee_out + file_name
                
                subprocess.call(["singularity", "exec", "-H", "/speed-scratch/shiva/",'/speed-scratch/bioinformatics-group/bioinformatics-singularity.beta.simg',"t_coffee", "-infile", tm_coffee_input_path, "-outfile", tm_coffee_out_put_path,
                                 "-evaluate", "-output", "fasta"])
                
        self.count_columns(tm_coffee_out)
        return tm_coffee_out

    def count_columns(self,column_files):
        d = os.walk(column_files, topdown=True)
        array=[]
        statics=[]

        for path, dir_list, file_list in d:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    array.append(file_name)
                else:
                    print('Invalid File Name :' + file_name)

            
            for filename in array:
                count=[]
                for f, seq_record in enumerate(SeqIO.parse(path+filename, 'fasta')):
                    tmp=0
                    for i in seq_record.seq:
                        if i !='-':
                            tmp+=1
                    count.append(tmp)
                if (len(count)!=0):
                    saved= str(int(sum(count)/len(count)+0.5))
                    org_length=str(len(seq_record.seq))
                    statics.append(': '.join([filename[:-6],saved+'/'+org_length ]) + '\n')
                

        self.writeToFile(statics, self.folder+self.folder[:-1]+".txt")

    def compute_xdet(self):
        xdet_output=self.detec_pos()
        # windows path
        # sds_xdet_result=str(Path().absolute())+"\SDS_XDet_result\"
        sds_xdet_result = self.folder+"SDS_XDet_result/"
        if os.path.exists(sds_xdet_result):
            shutil.rmtree(sds_xdet_result)
        os.mkdir(sds_xdet_result)
        self.startPoint(xdet_output,sds_xdet_result,'xdet')
        self.count_columns(sds_xdet_result)
        return sds_xdet_result

    def groupsim_preprocess(self):
        input_subclass = self.compute_secator()

        clean_input, subclasses = self.compute_subclass(input_subclass)

        groupsim_input=self.folder+"groupsim_input/"
        if os.path.exists(groupsim_input):
            shutil.rmtree(groupsim_input)
        os.mkdir(groupsim_input)
        c=os.walk(clean_input)
        list = []
        for path, dir_list, file_list in c:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    list.append(file_name)
                else:
                    print('Invalid File Name :' + file_name)
            for file_name in list:
                input_path = path + file_name
                filename = file_name[0:(len(file_name) - 6)]
                records = subclasses + filename
                out_put_path = groupsim_input + filename + '.clustal'
                fs = open(input_path, 'r+')
                gp = open(records, 'r')
                # lines = fs.readlines()
                groups = gp.readlines()
                count = 0
                tmp = ''
                line2 = fs.readline()
                for i in range(len(groups)):
                    for j in range(int(groups[i])):
                        if line2[:1] == ">":
                            line2 = line2.replace("\n", "")
                            count += 1
                            if count <= int(groups[i]):
                                tmp = line2[4:10] + "|group" + str(1 + i) + "    "

                                line2=fs.readline()
                        while(line2[:1] != ">" and line2!=''):
                            tmp = tmp + line2.replace("\n", "")
                            line2=fs.readline()
                        tmp = tmp+'\n'


                        f = open(out_put_path, 'a')
                        f.writelines(tmp)
                        # line2=tmp

                    count=0
        return groupsim_input,subclasses

    def run_groupsim(self):
        output = self.folder+"groupSim_result/"
        groupsim="group_sim_sdp.py"
        metric="blosum62.bla"
        if os.path.exists(output):
            shutil.rmtree(output)
        os.mkdir(output)

        input , subclasses= self.groupsim_preprocess()
        input=os.walk(input,topdown=True)
        list = []
        for path, dir_list, file_list in input:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    list.append(file_name)
                else:
                    print('Invalid File Name :' + file_name)
            for file_name in list:
                input_groupsim= path + file_name
                group_file_name = file_name[0:(len(file_name) - 8)]
                output_groupsim = output + group_file_name
                tmp = subclasses + group_file_name
                gp = open(tmp, 'r')
                records = len(gp.readlines())
                if records>1:

                    exe=["python2", groupsim, "-n", "-m", metric, "-w", "1", "-o",output_groupsim, input_groupsim]
                    for record in range(records):
                        group='group'+str(record+1)
                        exe.append(group)
                    subprocess.call(exe)
        return output

    def compute_groupsim(self):
        input = self.run_groupsim()
        groupsim_result = self.folder+"SDS_GroupSim_result/"
        if os.path.exists(groupsim_result):
            shutil.rmtree(groupsim_result)
        os.mkdir(groupsim_result)
        self.startPoint(input, groupsim_result, 'groupsim')
        self.count_columns(groupsim_result)
        return groupsim_result





    def run_speer(self):
        output=self.folder+"SpeerServer_result/"
        if os.path.exists(output):
            shutil.rmtree(output)
        os.mkdir(output)
        # speer_exe=str(Path().absolute())+'/tools/Speer/src/SPEER'
        input_subclass=self.compute_secator()
        clean_input,subclasses=self.compute_subclass(input_subclass)
        # clean_input='secator/clear_fasta/'
        # subclasses='secator/subclasses/'
        clean_input = os.walk(clean_input, topdown=True)
        list = []
        for path, dir_list, file_list in clean_input:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    list.append(file_name)
                else:
                    print('Invalid File Name :' + file_name)
            for file_name in list:
                input_speer = path + file_name
                speer_file_name = file_name[0:(len(file_name) - 6)]
                output_speer = output + speer_file_name 
                tmp = subclasses + speer_file_name

                gp = open(tmp, 'r')
                records = len(gp.readlines())
                if records > 1:
                    subprocess.call(["singularity", "exec", "-H", "/speed-scratch/shiva/",'/speed-scratch/bioinformatics-group/speer-new2.simg',"SPEER", "-r4sDir", "/usr/bin/","-wEDist", "1", "-wRE", "1", "-wERate", "1","-i", input_speer, "-o", output_speer, "-pf", tmp])
        return output


    def compute_speer(self):
        input_speer=self.run_speer()
        # input_speer="SpeerServer_result/"
        sds_speer_result = self.folder+"SDS_SpeerServer_result/"
        if os.path.exists(sds_speer_result):
            shutil.rmtree(sds_speer_result)

        os.mkdir(sds_speer_result)
        self.startPoint(input_speer,sds_speer_result,'speer')
        self.count_columns(sds_speer_result)
        return sds_speer_result



    def compute_secator(self):
        secator_out=self.folder+"secator/"
        secator_fasta = secator_out + "secator_fasta/"
        secator_clusters = secator_out + "secator_clu/"
        # secator_path = str(Path().absolute())+"/tools/Secator/secator"
        if os.path.exists(secator_out):
            shutil.rmtree(secator_out)
        os.mkdir(secator_out)
        os.mkdir(secator_clusters)
        os.mkdir(secator_fasta)
        list = []
        secator_input = os.walk(self.data, topdown=True)
        for path,dir_list,file_list in secator_input:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    list.append(file_name)
                else:
                    print ('Invalid File Name :' + file_name)
            for file_name in list:
                input_secator = path + file_name
                secator_file_name = file_name[0:(len(file_name)-6)]
                output_secator = secator_fasta + secator_file_name + '.fasta'
                output_cluster = secator_clusters + secator_file_name + '.clu'
                otfa="-otfa="+output_secator
                oclu="-oclu="+output_cluster
                subprocess.call(["singularity", "exec", "-H", "/speed-scratch/shiva/",'/speed-scratch/bioinformatics-group/bioinformatics-singularity.simg',"secator",input_secator,"-dt=alignment","-cm=hierar",otfa,oclu])

        return secator_fasta

    def compute_subclass(self,input):
        records = {}
        output=self.folder+'secator/subclasses/'
        if os.path.exists(output):
            shutil.rmtree(output)
        os.mkdir(output)
        output_input=self.folder+'secator/clear_fasta/'
        if os.path.exists(output_input):
            shutil.rmtree(output_input)
        os.mkdir(output_input)
        input=os.walk(input, topdown=True)
        list = []
        for path, dir_list, file_list in input:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    list.append(file_name)
                else:
                    print('Invalid File Name :' + file_name)
            for file_name in list:
                clear_file_name = file_name[0:(len(file_name)-6)]
                output_clear = output_input + clear_file_name + '.fasta'
                group = []
                tmp=output+str(clear_file_name)
                i=0
                seqs=[]
                for seq_record in SeqIO.parse(path+file_name, 'fasta'):
                    if (seq_record.id.startswith("GROUP_1")):
                        i = 0
                    else:
                        if(seq_record.id.startswith("GROUP")):
                            group.append(i)
                            i=0
                        else :
                            i+=1
                            seqs.append(seq_record)

                SeqIO.write(seqs, output_clear, "fasta")
                group.append(i)
                records[clear_file_name]=group

        for key in records.keys():
            with open(output+str(key), "a") as file:
                file.write("\n".join(str(e) for e in records[key]))
        return output_input,output


    def detec_pos(self):
        # windows path
        # metrix_path = str(Path().absolute())+"\metrix\blosum45.bla"
        # path_out = str(Path().absolute())+"\XDet_result\"
        metrix_path = "metrix/blosum45.bla"
        path_out = self.folder+"XDet_result/"
        if os.path.exists(path_out):
            shutil.rmtree(path_out)
        os.mkdir(path_out)
        c = os.walk(self.data, topdown=True)
        list = []
        for path,dir_list,file_list in c:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if result:
                    list.append(file_name)
                else:
                    print ('Invalid File Name :' + file_name)
            for file_name in list:
                input_path = path + file_name
                # filename = file_name[0:(len(file_name)-11)]
                out_put_path = path_out + file_name 
                with open(out_put_path, 'w') as outputFile:
                    # windows path
                    # xdet_exe=str(Path().absolute())+"\tools\JDet\programs\xdet_osx"
                    # xdet_exe = str(Path().absolute())+"/tools/JDet/programs/xdet_osx"
                    subprocess.call(["singularity", "exec", "-H", "/speed-scratch/shiva/",'/speed-scratch/bioinformatics-group/bioinformatics-singularity.beta.simg',"xdet", input_path, metrix_path, out_put_path, "-S", "10"], stdout=outputFile)
        return path_out


    def gp_readFileList(self,input):
        org_files = [each_file for each_file in os.listdir(self.folder+'secator/clear_fasta/')
                     if isfile(join(self.folder+'secator/clear_fasta/', each_file)) and not each_file.startswith('.')]

        out_files = [each_file for each_file in os.listdir(input)
                     if isfile(join(input, each_file)) and not each_file.startswith('.')]

        index_list = []
        if os.path.exists(self.folder+"groupsim_garbage/"):
            shutil.rmtree(self.folder+"groupsim_garbage/")
        os.mkdir(self.folder+"groupsim_garbage/")

        for i in range(len(org_files)):
            if org_files[i][:-6] not in out_files:
                index_list.append(join(self.folder+'secator/clear_fasta/', org_files[i]))
        for i in range(len(index_list)):
            shutil.move(index_list[i], self.folder+"groupsim_garbage/")

        org_files = [join(self.folder+'secator/clear_fasta/', each_file) for each_file in os.listdir(self.folder+'secator/clear_fasta/')
                     if isfile(join(self.folder+'secator/clear_fasta/', each_file)) and not each_file.startswith('.')]

        out_files = [join(input, each_file) for each_file in os.listdir(input)
                     if isfile(join(input, each_file)) and not each_file.startswith('.')]

        return org_files, out_files

    def speer_readFileList(self,input):
        org_files = [each_file for each_file in os.listdir(self.folder+'secator/clear_fasta/')
                     if isfile(join(self.folder+'secator/clear_fasta/', each_file)) and not each_file.startswith('.')]

        out_files = [each_file for each_file in os.listdir(input)
                     if isfile(join(input, each_file)) and not each_file.startswith('.')]

        index_list = []
        if os.path.exists(self.folder+"speer_garbage/"):
            shutil.rmtree(self.folder+"speer_garbage/")
        os.mkdir(self.folder+"speer_garbage/")

        for i in range(len(org_files)):
            if org_files[i][:-6] not in out_files:
                index_list.append(join(self.folder+'secator/clear_fasta/', org_files[i]))
        for i in range(len(index_list)):
            shutil.move(index_list[i], self.folder+"speer_garbage/")

        org_files = [join(self.folder+'secator/clear_fasta/', each_file) for each_file in os.listdir(self.folder+'secator/clear_fasta/')
                     if isfile(join(self.folder+'secator/clear_fasta/', each_file)) and not each_file.startswith('.')]

        out_files = [join(input, each_file) for each_file in os.listdir(input)
                     if isfile(join(input, each_file)) and not each_file.startswith('.')]

        return org_files, out_files

    def readFileList(self,input,type):
        if type=='groupsim':
            return self.gp_readFileList(input)
        elif type =='speer':
            return self.speer_readFileList(input)
        else:
            # get original input files for the speer server
            org_files = [join(self.data, each_file) for each_file in os.listdir(self.data)
                         if isfile(join(self.data, each_file)) and not each_file.startswith('.')]

            # get Xdet output files
            # windows path
            # out_files = [join(str(Path().absolute()) + "\XDet_result\", each_file) for each_file in
            #              os.listdir(str(Path().absolute()) + "\XDet_result\")
            #              if isfile(join(str(Path().absolute()) + "\XDet_result\", each_file)) and not each_file.startswith(
            #         '.')]
            out_files = [join(input, each_file) for each_file in os.listdir(input)
                         if isfile(join(input, each_file)) and not each_file.startswith('.')]



            if len(org_files) != len(out_files):
                # print("The number of input files is not same as the number of output files")
                sys.exit(0)

            return org_files, out_files

    def readOriInput(self,filename):
        recordList = []
        temp = []
        key = ''
        with open(filename, 'r') as inputfile:
            oneline = inputfile.readline()
            while oneline:
                if oneline.startswith('>'):
                    if not temp:
                        key = oneline.replace('\n', '')
                    else:
                        recordList.append({key: ''.join(temp)})
                        key = oneline.replace('\n', '')
                        temp = []
                else:
                    temp.append(oneline.replace('\n', ''))
                oneline = inputfile.readline()

        recordList.append({key: ''.join(temp)})
        return recordList

    def speer_readOutputFile(self,filename):

        positions = []

        with open(filename, 'r') as outputfile:
            oneline = outputfile.readline()
            while oneline:
                parts = oneline.replace('\n', '').split()
                if len(parts) == 5 and parts[0].isdigit():
                    positions.append(int(parts[0]))
                oneline = outputfile.readline()

        positions = set(positions)
        return positions

    def xdet_readOutputFile(self,filename):
        positions = []
        
        with open(filename, 'r') as outputfile:
            oneline = outputfile.readline()
            if oneline =='\n':
                return set([])
            
            while oneline:
                parts = oneline.replace('\n', '').split()
                
                positions.append(int(parts[0]))
                oneline = outputfile.readline()

        positions = set(positions)
        return positions

    def goupsim_readOutputFile(self,filename):
        positions = []

        with open(filename, 'r') as outputfile:
            oneline = outputfile.readline()
            while oneline:
                parts = oneline.replace('\n', '').split()
                if len(parts) == 5 and parts[0].isdigit() and parts[1] != 'None':
                    positions.append(int(parts[0]))
                oneline = outputfile.readline()

        positions = set(positions)
        return positions


    def readOutputFile(self,filename,type):
        if type=='speer':
            return self.speer_readOutputFile(filename)
        elif type=='xdet':
            return self.xdet_readOutputFile(filename)
        elif type=='groupsim':
            return self.goupsim_readOutputFile(filename)

    def writeToFile(self,output, filename):
        # windows path
        # filename = str(Path().absolute()) + "\SDS_XDet_result\" + filename

        with open(filename, 'w') as outputfile:
            for eachData in output:
                outputfile.writelines(eachData)
            outputfile.flush()

    def startPoint(self,input,output_folder,type):
        path_splitor = "/"
        output = []
        dashStat = []
        original_files, output_files = self.readFileList(input,type)
        for i in range(len(original_files)):
            positions = self.readOutputFile(output_files[i],type)
            recordList = self.readOriInput(original_files[i])
            counter = 0
            length = 0
            counting = True
            for eachRecord in recordList:
                original = ['dummy']
                item = list(eachRecord)[0]
                key = item[0:len(item)]
                value = list(eachRecord.values())[0]
                original += list(value)
                target = []

                for j in range(len(original)):
                    if j in positions:
                        target.append(original[j])
                        if counting:
                            counter += 1
                    else:
                        target.append('-')

                output.append({key: ''.join(target)})

                if counting:
                    length = (len(original) - 1)
                    counting = False

            toFile = []
            for element in output:
                item = list(element)[0]
                toFile.append('\n'.join([item[0:len(item)], list(element.values())[0], '']))

            self.writeToFile(toFile, output_folder+original_files[i].split(path_splitor)[-1])
            dashStat.append({original_files[i].split(path_splitor)[-1]: '/'.join([str(counter), str(length)])})

            # print("Replace characters to '-' is done for %s" % original_files[i].split(path_splitor)[-1])
            output = []

        statistics = []
        for element in dashStat:
            item = list(element)[0]
            statistics.append(': '.join([item[0], list(element.values())[0]]) + '\n')
        # windows path
        # self.writeToFile(statistics, str(Path().absolute())+"\SDS_XDet_result\dash_statistics.txt")
        self.writeToFile(statistics, output_folder+"/dash_statistics.txt")