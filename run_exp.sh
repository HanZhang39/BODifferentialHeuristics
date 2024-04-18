
mkdir den520d
mkdir den520d/0
mkdir den520d/001

# precomputed landmarks are in the assets/landmarks folder. You can check ./scripts/compute_landmarks.ipynb for the scrip used to compute them.

while read p; do
    echo $p
    # compute single to all Pareto frontiers
    ./bin/bod --start $p --map assets/maps/den520d_random_10_0.gr assets/maps/den520d_random_10_1.gr -o den520d/0/den520d_output_${p}.txt -a BOD


    # compress single to all Pareto frontiers with an epsilon-value of 0.01
    ./bin/compress -i den520d/0/den520d_output_${p}.txt -e 0.01 -o den520d/001/den520d_${p}_apex.txt -p den520d/001/den520d_${p}_path.txt
done < assets/landmarks/landmarks_den520d.txt


# solver with differential heuristics
./bin/solver_dh -q ./assets/instances/den520d_query.txt -m assets/maps/den520d_random_10_0.gr assets/maps/den520d_random_10_1.gr -o ./output_dh.txt --dh den520d_001_apex.txt --dhp den520d_001_path.txt --nl 128 --update_interval 500 --update_threshold 0.05

# baseline solver 
./bin/solver -q ./assets/instances/den520d_query.txt -m assets/maps/den520d_random_10_0.gr assets/maps/den520d_random_10_1.gr -o ./output.txt -a NAMOA
