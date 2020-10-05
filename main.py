from final_pipeline import FinalPipeline
import sys


def main():
    train = str(sys.argv[1])
    test = str(sys.argv[2])
    database = str(sys.argv[3])
    msa = str(sys.argv[4]).lower()
    sds = str(sys.argv[5]).lower()
    output = FinalPipeline(train, test, database, msa, sds)
    output.result()
    print("done!")


if __name__=="__main__":
    main()
