using Munkres, Combinatorics, Clustering

function evaluate_partition(true_blocks::Vector{Int64}, alg_blocks::Vector{Int64})
    true_blocks_set = Set(true_blocks)
    alg_blocks_set = Set(alg_blocks)
    if -1 in true_blocks_set
        pop!(true_blocks_set, -1) # -1 is unknown
    end
    println("")
    println("Partition Correctness Evaluation")
    println("Number of nodes: $(length(alg_blocks))")
    println("Number of partitions in truth partition: $(length(true_blocks_set))")
    println("Number of partitions in alg. partition: $(length(alg_blocks_set))")
    println("")

    confusion_matrix = zeros(length(true_blocks_set), length(alg_blocks_set))
    for i=1:length(alg_blocks)
        if true_blocks[i] != -1
            confusion_matrix[true_blocks[i], alg_blocks[i]] += 1
        end
    end

    accuracy = compute_accuracy!(confusion_matrix)

    pairwise_precision, pairwise_recall = pairwise_metrics(true_blocks, alg_blocks)
    adj_rand_index, rand_index, _, _  = randindex(true_blocks, alg_blocks)

    println("Rand Index: $rand_index")
    println("Adjusted Rand Index: $adj_rand_index")
    println("")

    accuracy, pairwise_precision, pairwise_recall, adj_rand_index, rand_index, information_theory_metrics(confusion_matrix)
end


function compute_accuracy!(confusion_matrix::Array{Float64, 2})
    if size(confusion_matrix, 1) > size(confusion_matrix, 2)
        confusion_matrix = transpose(confusion_matrix)
    end

    best_assignments = munkres(-confusion_matrix)

    confusion_matrix_before_assignment = deepcopy(confusion_matrix)
    total = 0
    for (row, col) in enumerate(best_assignments)
        confusion_matrix[:, row] = confusion_matrix_before_assignment[:, col]
        total += confusion_matrix[row, row]
    end

    unassociated_col = setdiff(Set(1:size(confusion_matrix, 2)), Set(best_assignments))
    counter = 0
    for col in unassociated_col
        confusion_matrix[:, size(confusion_matrix, 1) + counter] = confusion_matrix_before_assignment[:, col]
        counter += 1
    end

    if size(confusion_matrix, 2) > size(confusion_matrix, 1)
        confusion_matrix = transpose(confusion_matrix)
    end

    # joint probability of the two partitions is just the normalized contingency table
    joint_prob = confusion_matrix / sum(confusion_matrix)
    accuracy = sum(diag(joint_prob))

    println("Contingency Table: ", confusion_matrix)
    println("Accuracy (with optimal partition matching): $accuracy")

    accuracy
end

function pairwise_metrics(true_blocks::Vector{Int64}, alg_blocks::Vector{Int64})
    counts = zeros(Int64, 4)
    for (n1, n2) in combinations(1:length(alg_blocks), 2)
        # Same truth block
        if true_blocks[n1] == true_blocks[n2]
            if alg_blocks[n1] == alg_blocks[n2]
                counts[1] += 1
            else
                counts[3] += 1
            end
        # Different truth blocks
        else
            if alg_blocks[n1] == alg_blocks[n2]
                counts[4] += 1
            else
                counts[2] += 1
            end
        end
    end
    pairwise_precision = counts[1]/(counts[1] + counts[4])
    pairwise_recall = counts[1]/(counts[1] + counts[3])
    println("Pairwise Recall: $pairwise_recall")
    println("Pairwise Precision: $pairwise_precision")
    pairwise_precision, pairwise_recall
end

function information_theory_metrics(confusion_matrix::Array{Float64, 2})
    joint_prob = confusion_matrix / sum(confusion_matrix)
    marginal_prob_b2 = sum(joint_prob, 1)
    marginal_prob_b1 = sum(joint_prob, 2)
    idx1, _ = findn(marginal_prob_b1)
    _, idx2 = findn(marginal_prob_b2)
    conditional_prob_b2_b1 = zeros(size(joint_prob))
    conditional_prob_b1_b2 = zeros(size(joint_prob))
    conditional_prob_b2_b1[idx1, :] = joint_prob[idx1, :] ./ marginal_prob_b1[idx1]
    conditional_prob_b1_b2[:,idx2] = joint_prob[:, idx2] ./ marginal_prob_b2[idx2]
    H_b2 = -sum(marginal_prob_b2[idx2] .* log.(marginal_prob_b2[idx2]))
    H_b1 = -sum(marginal_prob_b1[idx1] .* log.(marginal_prob_b1[idx1]))
    rows, cols = findn(joint_prob)
    H_b2_b1 = -sum([joint_prob[r, c] * log(conditional_prob_b2_b1[r, c]) for (r, c) in  zip(rows, cols)])
    H_b1_b2 = -sum([joint_prob[r, c] * log(conditional_prob_b1_b2[r, c]) for (r, c) in  zip(rows, cols)])

    marginal_prod = marginal_prob_b1 * marginal_prob_b2

    MI_b1_b2 = sum([joint_prob[r, c] * log(joint_prob[r, c] / marginal_prod[r, c]) for (r, c) in  zip(rows, cols)])

    if H_b1>0
        fraction_missed_info = H_b1_b2 / H_b1
    else
        fraction_missed_info = 0
    end
    if H_b2>0
        fraction_err_info = H_b2_b1 / H_b2
    else
        fraction_err_info = 0
    end

    println("Entropy of truth partition: ", abs(H_b1))
    println("Entropy of alg. partition: ", abs(H_b2))
    println("Conditional entropy of truth partition given alg. partition: ", abs(H_b1_b2))
    println("Conditional entropy of alg. partition given truth partition: ", abs(H_b2_b1))
    println("Mututal information between truth partition and alg. partition: ", abs(MI_b1_b2))
    println("Fraction of missed information: ", abs(fraction_missed_info))
    println("Fraction of erroneous information: ", abs(fraction_err_info))

    H_b1, H_b2, H_b1_b2, H_b2_b1, MI_b1_b2, fraction_missed_info, fraction_err_info
end
