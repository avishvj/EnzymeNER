import argparse

import tabolistem_model as tabm

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='PROG')


    parser.add_argument('-t', '--text_path', type=str)
    parser.add_argument('-a', '--annot_path', type=str)
    parser.add_argument('-o', '--output_name', type=str)

    args = parser.parse_args()

    tm = tabm.TaboListem()
    #tm.load('tabolistem_TaboListemModel.json','./TrainedModels/epoch_9_TaboListemModel_weights')
    #result = tm.process('ELISAs were performed with 96-well microtiter plates (Corning Costar, Cambridge, Mass.) coated with Mtb81, 38-kDa antigen (200 ng/well), or TB lysate (100 ng/well). Coating was done overnight at 4\u00b0C. Plates were then aspirated and blocked with PBS containing 1% (wt/vol) bovine serum albumin for 2 h at room temperature, followed by a wash with PBS containing 0.1% Tween 20 (PBST2). Serum (diluted 1/25 in PBST2 for Mtb81 and 1/100 for the 38-kDa antigen and the lysate) was added to the wells and incubated for 30 min at room temperature. Following incubation, the wells were washed six times with PBST2 and incubated with protein A-horseradish peroxidase conjugate (Sigma Chemical Co., St. Louis, Mo.) at a 1/20,000 dilution for 30 min. Plates were then washed six times with PBST2 and incubated with tetramethylbenzidine substrate (Kirkegaard & Perry Laboratories, Inc., Gaithersburg, Md.) for a further 15 min. The reaction was stopped by the addition of 1 N sulfuric acid, and plates were read at 450 nm with an ELISA plate reader (Biotek, Hyland Park, Va.). The cutoff for the assays was the mean of the negative population plus three standard deviations of the mean.')

    #print(result)
    #tm.train('TestingSet.txt', "TestingSetAnnot.txt", "TaboListemModel")
    tm.load('tabolistem_BioBERT_TaboListemModel.json','./BioBertModels/epoch_8_TaboListemModel_weights')
    #result = tm.process('Extracellular-exposed caveolae-specific proteins CD36 and copper-containing amine oxidase were concealed inside the vesicles and resisted trypsin treatment.')
    #'Kinase is a kind of enzyme.')
    #print(result)
    tm.test('GoldSet.txt', "GoldSetAnnot.txt")
    #tm.train(args.text_path, args.annot_path,
     #        args.output_name)

# e.g. python run_app.py -t "TrainingSet.txt" -a "TrainingSetAnnot.tsv" -o "TaboListemModel"c
# e.g. python run_app.py -t "TrainingSet.txt" -a "TrainingSetAnnot.tsv" -o "TaboListemModel"c
