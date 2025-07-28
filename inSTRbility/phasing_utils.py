import numpy as np


def build_cluster_haplotypes(cluster, significant_snps, SNPS):
    cluster_haplotype = {}
    for read_idx in cluster:
        for snp_pos in significant_snps:
            alleles = list(filter(lambda x: x not in {'COV', 'Q'}, SNPS[snp_pos].keys()))
            for allele in alleles:
                if read_idx in SNPS[snp_pos][allele]:
                    if snp_pos not in cluster_haplotype:
                        cluster_haplotype[snp_pos] = {allele: 1}
                    else:
                        if allele not in cluster_haplotype[snp_pos]:
                            cluster_haplotype[snp_pos][allele] = 1
                        else:
                            cluster_haplotype[snp_pos][allele] += 1
    for snp_pos in cluster_haplotype:
        cluster_haplotype[snp_pos] = sorted(cluster_haplotype[snp_pos].items(), key=lambda item: item[1], reverse=True)[0][0]
    return cluster_haplotype

def print_cluster_haplotypes(clusters, read_indices, READ_NAMES):
    snps = []
    for cluster in clusters:
        snps +=  list(cluster.keys())
    snps = sorted(set(snps))
    for i, cluster in enumerate(clusters):
        alleles = []
        for snp in snps:
            if snp in cluster:
                alleles.append(cluster[snp])
            else:
                alleles.append('-')
        print(READ_NAMES[read_indices[i]], ''.join(alleles), sep='\t')

def phase_unphased_reads(unphased_read, snps, final_haplotypes, SNPS):
    """
    phase the unphased reads based on least changing Q value
    Args:
        unphased_reads (set): Set of read indices that are not phased.
        final_haplotypes (list): List of tuples containing the phased reads.
        filtered_snp_info (dict): Dictionary containing SNP positions and their alleles with supporting reads.
        SNPS (dict): Dictionary containing SNP information.
        READ_VARS (dict): Dictionary containing read variables.
        LOCI_VARS (dict): Dictionary containing locus variables.
    Returns:
        None: Updates the final haplotypes with unphased reads.
    """

    # build snp haplotypes for the reads withing the phased groups
    # for each uphased read we check the alleles at the SNP positions and assess sum of the Qvalues where the read is present

    cluster1, cluster2 = final_haplotypes
    cluster1_haplotype = build_cluster_haplotypes(cluster1, snps, SNPS)
    cluster2_haplotype = build_cluster_haplotypes(cluster2, snps, SNPS)

    phased_reads = set()

    for read_idx in unphased_read:
        # read_snps = set(READ_VARS[read_idx]['X'])
        read_haplotype = {}
        for snp_pos in snps:
            alleles = list(filter(lambda x: x not in {'COV', 'Q'}, SNPS[snp_pos].keys()))
            for allele in alleles:
                if read_idx in SNPS[snp_pos][allele]:
                    if allele != 'r':
                        read_haplotype[snp_pos] = {allele: SNPS[snp_pos]['Q'][read_idx]}
                    else: read_haplotype[snp_pos] = { allele: 30 }

        c1_mismatch = 0; c2_mismatch = 0
        c1_match = 0; c2_match = 0
        informative_snps = 0
        for snp_pos in read_haplotype:
            read_allele = list(read_haplotype[snp_pos].keys())[0]
            c1_allele = cluster1_haplotype[snp_pos] if snp_pos in cluster1_haplotype else None
            c2_allele = cluster2_haplotype[snp_pos] if snp_pos in cluster2_haplotype else None
            if c1_allele is not None and read_allele != c1_allele:
                c1_mismatch -= read_haplotype[snp_pos][read_allele]
                informative_snps += 1
            if c2_allele is not None and read_allele != c2_allele:
                c2_mismatch -= read_haplotype[snp_pos][read_allele]
                informative_snps += 1
            if c1_allele is not None and read_allele == c1_allele:
                c1_match += 1
            if c2_allele is not None and read_allele == c2_allele:
                c2_match += 1

        phased = False
        if informative_snps > 0:
            if   c1_mismatch > c2_mismatch and c1_match > 0:
                phased = True; final_haplotypes[0].add(read_idx)
            elif c2_mismatch > c1_mismatch and c2_match > 0:
                phased = True; final_haplotypes[1].add(read_idx)

        if phased: phased_reads.add(read_idx)

    for read_id in phased_reads:
        unphased_read.remove(read_id)


def delete_lowcov_alleles(passed_snps, filtered_snp_info, args):

    # for each snp position we get rid of the nucs that have coverage lesser than user parameters
    for position in passed_snps:
        position_cov = sum([len(reads) for reads in filtered_snp_info[position]['alleles'].values()]) # calculate the total reads for that snp

        del_alleles = []
        for allele in filtered_snp_info[position]['alleles']:
            fraction_allele_reads = len(filtered_snp_info[position]['alleles'][allele]) / position_cov
            if fraction_allele_reads < args.snp_reads: del_alleles.append(allele) # if the reads in 'nuc' have reads less than 25%, delete it

        for allele in del_alleles:
            del filtered_snp_info[position]['alleles'][allele]

    # removing SNP positions that have less than 2 nucleotides alleles after filtering
    del_positions = []
    for pos in filtered_snp_info:
        if len(filtered_snp_info[pos]['alleles']) < 2:
            del_positions.append(pos)
    for pos in del_positions:
        del filtered_snp_info[pos]


def haplocluster_reads(sorted_snps, snp_info, read_indices, SNPS, READ_NAMES, args):
    """
    Function to cluster the reads based on SNPs and their supporting reads.

    Args:
        snp_allelereads (dict): Dictionary containing SNP positions and their alleles with supporting reads.
        coverage_sorted_snps (list): List of SNP positions sorted by coverage.
        read_indices (set): Set of read indices that support the SNPs.
        args (Namespace): Command line arguments containing parameters like snpC, snpR, phasingR, snpQ.

    Returns:
        list: A list containing the final haplotypes, minimum SNP position, skip point,
    """

    threshold_range = [(0.3, 0.7),(0.25, 0.75),(0.2, 0.8)] # threshold values to get Significant_snps
    for idx, range in enumerate(threshold_range):

        final_haplotypes = ()
        fail_code = 10

        r1 = range[0] # threshold value 1
        r2 = range[1] # threshold value 2

        filtered_snp_info = {}
        passed_snps     = []

        for pos in sorted_snps:
            # if the position is covered by less than 60% of the reads corresponding to the locus, skip it
            if snp_info[pos]['cov'] < (0.6*len(read_indices)):
                fail_code = 3; break

            passed_alleles = 0
            pos_coverage = snp_info[pos]['cov']
            for allele in snp_info[pos]['alleles']:
                allele_cov = len(snp_info[pos]['alleles'][allele])
                # allele coverage falls within the range for possible segregation
                if r1 * pos_coverage <= allele_cov <= r2 * pos_coverage:
                    passed_alleles += 1

            if passed_alleles >= 2: # if the SNP position has atleast 2 nucs that pass the threshold
                passed_snps.append(pos)
                filtered_snp_info[pos] = {'cov': snp_info[pos]['cov'], 'alleles': snp_info[pos]['alleles'], 'Q': snp_info[pos]['Q']}

        if len(filtered_snp_info) == 0:
            if idx < 2:   # move to the next threshold range
                continue
            else:
                fail_code = 2

                return [final_haplotypes, fail_code, 0]

        if len(passed_snps) == 0:
            fail_code = 2
            return [(), fail_code, len(passed_snps)]

        passed_snps = passed_snps[:args.snp_count] # limiting the number of SNPs to be considered for phasing

        delete_lowcov_alleles(passed_snps, filtered_snp_info, args) # remove alleles which have low coverage in the passed SNPs

        # if all the SNP positions failed the coverage filtering then we skip phasing
        if len(filtered_snp_info) == 0: # after removing snp based on their nuc's read coverage, if there are no snp left then go to next threshold range
            fail_code = 2
            return [(), fail_code, len(passed_snps)]

        result = merge_snpreadsets(passed_snps, filtered_snp_info, read_indices, sorted_snps, SNPS, READ_NAMES, args) # calling the phasing function
        final_haplotypes, status, fail_code, snp_num = result

        if status: break

        if idx == 2: #level_split:
            break

    return [final_haplotypes, fail_code, snp_num]


def merge_snpreadsets(snps, snp_info, read_indices, unfiltered_snps, SNPS, READ_NAMES, args):
    """
    Merge the reads sets to create two final haplotype groups based on SNPs
    Args:
        filtered_snps (dict): Dictionary containing SNP positions and their alleles with supporting reads.
        ordered_split_pos (list): List of SNP positions sorted by coverage.
        read_indices (set): Set of read indices that support the SNPs.
        args (Namespace): Command line arguments containing parameters like snpC, snpR, phasingR, snpQ.
    Returns:
        tuple: A tuple containing the final haplotypes, status (True/False), minimum SNP position, skip point,
               chosen SNP quality values, phased reads, and the number of SNP
    """

    fail_code = 10

    # SNP Combination method
    snps = list(filter(lambda x: x in snp_info, snps)) # filter the snps to get only those which are present in snp_info

    read_clusters = {}
    for idx in range(len(snps)): # Comparison starting from 1st snp to all other snp

        # if we have reached the last SNP position and there is no read cluster, then we break
        if (read_clusters != {}) and idx == len(snps)-1: break

        snp_A = snps[idx] # current SNP position
        readsets_A = list(snp_info[snp_A]['alleles'].values()) # [ {1,2,3,4,5},  {6,7,8,9,10} ]  nuc's read set for current Target SNP

        # holds the number of reads common between the two SNP haplogroups
        read_clusters[snp_A] = {}
        for snp_B in snps[idx + 1:]:
            uniqueness = 0

            readsets_B = snp_info[snp_B]['alleles'].values()
            for set_B in readsets_B:  # [  {1,2,3,4,6},  {5,7,8,9,10}  ]  nuc's read set for Query SNP
                common_reads = set()
                for set_A in readsets_A:
                    intersection = set_A & set_B  # intersection is calculated for each combination of 4 read sets
                    common_reads.add(len(intersection))   # taking the min value from a pairwise interection values

                # storing the least number of common reads between the nucleotide read sets for the two SNP positions
                uniqueness += min(common_reads)

            read_clusters[snp_A][snp_B] = uniqueness # pos_cluster = { 1023 : {1036 : 0, 1045 : 1, 1123 : 3,..},  1036 : { 1045: 1, 1123 : 2,..}, ..............}

    significant_snps = []
    for snp_A in read_clusters:
        # this checks for exclusivity between the read groups of different nucleotide alleles at the SNP position
        if list(read_clusters[snp_A].values()).count(0) >= 2: # check whether any of the pos_cluster have atleast 2 zeros in their match scores
            significant_snps.append(snp_A)
            significant_snps.extend(sorted(read_clusters[snp_A].keys(), key=lambda item: read_clusters[snp_A][item])) # if yes take that pos_cluster and proceed for clustering
            break

    # if there are no significant SNPs, then we take the SNP position with least mismatch score
    if significant_snps == []: # if there are no pos_cluster with atleast 2 zeros in it
        most_unique = {}
        for snp_A in read_clusters:
            most_unique[snp_A] = sum(sorted(val for val in read_clusters[snp_A].values())[:2]) # take sum of 1st 2 mismatch scores from sorted pos:mis_score and append in to 'least_mismatches' dictionary
        most_unique_snp = sorted(most_unique.keys(), key=lambda item: most_unique[item])[0] # now take the snp position with least score and proceed for clustering
        significant_snps.append(most_unique_snp)
        significant_snps.extend(sorted(read_clusters[most_unique_snp].keys(), key=lambda item: read_clusters[most_unique_snp][item])) # [1023, 1036, 1045, 1123]

    significant_snps = list(filter(lambda x: x in snp_info, significant_snps)) # filtering the significant_snps to get only those which are present in final_ordered_dict

    

    cluster1, cluster2 = set(), set()
    for snp in significant_snps: # clustering starts from the position of how its arranged in the dict, The 2 read_set will be two new haplotypes
        for allele_reads in snp_info[snp]['alleles'].values():
            if not cluster1:
                cluster1 |= allele_reads
            elif not cluster2:
                cluster2 |= allele_reads

            elif (len(allele_reads & cluster1)>0.7*len(allele_reads)) and (len(allele_reads & cluster2)<0.05*len(allele_reads)): # then successive read_set are joined based on their intersection percentage
                cluster1 |= allele_reads
            elif (len(allele_reads & cluster2)>0.7*len(allele_reads)) and (len(allele_reads & cluster1)<0.05*len(allele_reads)):
                cluster2 |= allele_reads

    phased_reads   = cluster1 | cluster2
    unphased_reads = read_indices - phased_reads

    # if len(unphased_reads) > 0:
    #     print(unphased_reads)
    #     phase_unphased_reads(unphased_reads, unfiltered_snps, (cluster1, cluster2), SNPS)

    # after clustering the total reads in the both reads should be greater thn 50% of total supporting reads for that specific locus
    if ((len(cluster1)+len(cluster2)) >= args.fraction_phased_reads*len(read_indices)):
        # print('Cluster 1:', sorted(cluster1))
        # print('Cluster 2:', sorted(cluster2))
        # qvalue_phasing(significant_snps, read_indices, SNPS) # calling the qvalue phasing function to cluster the reads based on SNPs and their supporting reads

        return [(cluster1, cluster2), True, fail_code, len(snps)]
    else:
        fail_code = 4
        return [(), False, fail_code, len(snps)]


def score_cluster_membership(cluster, read_idx, read_indices, score_matrix):
    """
    Function to score the membership of a read in a cluster based on its supporting reads.

    Args:
        cluster (set): Set of read indices in the cluster.
        read_idx (int): Index of the read to be scored.
        read_indices (set): Set of all read indices.
        score_matrix (list): List of scores for each read pair.

    Returns:
        int: Score indicating the membership of the read in the cluster.
    """
    scores = []
    for read in cluster:
        i = read_indices.index(read)
        j = read_indices.index(read_idx)
        index = (i * len(read_indices)) + j
        scores.append(score_matrix[index])
    return np.mean(scores) if scores else 0



def qvalue_phasing(sorted_snps, read_indices, SNPS, haplotypes=None):
    """
    Function to cluster the reads based on SNPs and their supporting reads.

    Args:
        snp_allelereads (dict): Dictionary containing SNP positions and their alleles with supporting reads.
        coverage_sorted_snps (list): List of SNP positions sorted by coverage.
        read_indices (set): Set of read indices that support the SNPs.
        args (Namespace): Command line arguments containing parameters like snpC, snpR, phasingR, snpQ.

    Returns:
        list: A list containing the final haplotypes, minimum SNP position, skip point,
    """

    sorted_reads = sorted(read_indices)

    clusters = []
    for read_idx in read_indices:
        clusters.append(build_cluster_haplotypes({read_idx}, sorted_snps, SNPS))

    score_matrix = []
    for i, read_a in enumerate(sorted_reads):
        for j, read_b in enumerate(sorted_reads):
            if read_a == read_b: score_matrix.append(0); continue
            if j < i: score_matrix.append(score_matrix[j*len(sorted_reads) + i]); continue
            score = 0
            for snp_pos in sorted_snps:
                allele_a = None; qual_a = None
                allele_b = None; qual_b = None
                for allele in SNPS[snp_pos]:
                    if allele == 'COV' or allele == 'Q': continue
                    if read_a in SNPS[snp_pos][allele]:
                        allele_a = allele
                        if allele != 'r': qual_a = SNPS[snp_pos]['Q'][read_a]
                        else: qual_a = 25
                    if read_b in SNPS[snp_pos][allele]:
                        allele_b = allele
                        if allele != 'r': qual_b = SNPS[snp_pos]['Q'][read_b]
                        else: qual_b = 25
                if allele_a is None or allele_b is None:
                    pass
                elif allele_a == allele_b:
                    pass
                else:
                    score -= (qual_a + qual_b)
            score_matrix.append(score)

    sorted_indices = np.argsort(score_matrix)

    num_reads = len(read_indices)
    cluster1 = set()
    cluster2 = set()
    phased_reads = set()
    for itr, idx in enumerate(sorted_indices):
        if score_matrix[idx] == 0: break
        read_a = sorted_reads[idx // num_reads]
        read_b = sorted_reads[idx % num_reads]

        if read_a in phased_reads and read_b in phased_reads: continue

        if itr == 0:
            cluster1.add(read_a)
            cluster2.add(read_b)
            phased_reads.add(read_a); phased_reads.add(read_b)
            continue
        if read_b not in phased_reads:
            score_b1 = score_cluster_membership(cluster1, read_b, sorted_reads, score_matrix)
            score_b2 = score_cluster_membership(cluster2, read_b, sorted_reads, score_matrix)
            if   score_b1 < score_b2:
                cluster2.add(read_b)
            elif score_b1 > score_b2:
                cluster1.add(read_b)
            phased_reads.add(read_b)
        if read_a not in phased_reads:
            score_a1 = score_cluster_membership(cluster1, read_a, sorted_reads, score_matrix)
            score_a2 = score_cluster_membership(cluster2, read_a, sorted_reads, score_matrix)
            if   score_a1 < score_a2:
                cluster2.add(read_a)
            elif score_a1 > score_a2:
                cluster1.add(read_a)
            phased_reads.add(read_a)

        if len(phased_reads) == len(read_indices):
            break

    for read_a in read_indices:
        if read_a not in phased_reads:
            score_a1 = score_cluster_membership(cluster1, read_a, sorted_reads, score_matrix)
            score_a2 = score_cluster_membership(cluster2, read_a, sorted_reads, score_matrix)
            if   score_a1 < score_a2: cluster2.add(read_a)
            elif score_a1 > score_a2: cluster1.add(read_a)
            else: print('Unassigned read:', read_a)
            phased_reads.add(read_a)

    print('Cluster 1:', sorted(cluster1))
    print('Cluster 2:', sorted(cluster2))
