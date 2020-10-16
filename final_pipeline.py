from division import DIVISION
from making_database import DATABASE
from blastrun import BlastRun
from msa import MSA
from sds import SDS
from hmm import HMM
from latex import ReadResult
import os
import shutil

class FinalPipeline:

    def __init__(self,train,test,database,msa,sds):
        self.train=train
        self.test=test
        self.database=database
        self.msa=msa
        self.sds=sds
       



    def result(self):
        folder=self.msa+'-'+self.sds+'/'
        if os.path.exists(folder):
            shutil.rmtree(folder)
        os.makedirs(folder)

        
        training_set = DIVISION(self.train, "train",folder)
        training_set.divide()
        testing_set = DIVISION(self.test, "test",folder)
        testing_set.divide()
        database = DATABASE(self.train,self.test, self.database,folder)
        database.compute()
        blastrun = BlastRun(folder)
        msa_input = blastrun.compute()
        msa = MSA(self.msa, msa_input, folder)
        sds_input = msa.compute()
        sds = SDS(self.sds, sds_input,folder)
        hmm_input = sds.compute()
        hmm = HMM(hmm_input,folder)
        hmm.compute()
        final_result = ReadResult(folder)
        final_result.read_result()



