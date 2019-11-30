import os
import subprocess

executable = 'p2pRaster.out'
currDir = os.getcwd()
outLogFile = currDir + '/log.csv'
outDir = currDir + '/log' #Output directory for algorithms' output
if not os.path.exists(outDir):
    os.makedirs(outDir)

subprocess.call(["make"], stdout=subprocess.PIPE, cwd=currDir)

# inFiles = ['Iris.csv', 'HTRU_2']
# graphType = [2, 3] # 2 is Barabasi-Albert, 3 is Erdos-Renyi
# peers = [1000, 5000, 10000, 15000, 20000]
# default_peers = 10000
# fanOut = [1, 3, -1]
# default_fanOut = 1
# rounds = [20, 22, 24, 26, 28]
# default_rounds = 24

inFiles = ['s1.csv']
graphType = [2] # 2 is Barabasi-Albert, 3 is Erdos-Renyi
peers = [10]
fanOut = [3]
rounds = [5]

totRuns = 0
for input in inFiles:
    datasetPath = currDir + '/Datasets/S-sets/' + input
    for graph in graphType:
        for npeers in peers:
            totRuns += 1
            outFile = input + '-d_' + str(graph) + '-p_' + str(npeers)
            subprocess.call([currDir + '/' + executable,
                                   "-p", str(npeers),
                                   "-d", str(graph),
                                   "-i", datasetPath,
                                   "-o", outFile])
        for fOut in fanOut:
            totRuns += 1
            outFile = input + '-d_' + str(graph) + '-f_' + str(fOut)
            subprocess.call([currDir + '/' + executable,
                                   "-f", str(fOut),
                                   "-d", str(graph),
                                   "-i", datasetPath,
                                   "-o", outFile])
        for nround in rounds:
            totRuns += 1
            outFile = input + '-d_' + str(graph) + '-n_' + str(nround)
            subprocess.call([currDir + '/' + executable,
                                   "-r", str(nround),
                                   "-d", str(graph),
                                   "-i", datasetPath,
                                   "-o", outFile])
print(totRuns)

#subprocess.call(["make", "clean"])
