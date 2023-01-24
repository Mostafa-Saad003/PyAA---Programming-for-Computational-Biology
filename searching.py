def Naive(pattern, sequence):
    pattern_length = len(pattern)  # corresponds to window size
    sequence_length = len(sequence)
    list_of_matches = []
    # loops over sequence
    for i in range(sequence_length-pattern_length+1):
        j = 0
        # loops over pattern and sequence simultaneously to check for matches between the pattern and the window in question
        while j < pattern_length:
            if sequence[i+j] != pattern[j]:
                break
            j += 1
            # if a perfect match is found between the pattern and the full window, the indices corresponding to the matching window are reported
            if j == pattern_length:
                list_of_matches.append((i+1,i+j))
    return list_of_matches

def BoyerMoore(pattern, sequence):
    # Preprocessing
    m = len(pattern)
    n = len(sequence)
    bad_char = [-1] * 256
    patterns = []
    for i in range(m):
        bad_char[ord(pattern[i])] = i
    # Searching
    i = 0
    while i <= n - m:
        j = m - 1
        while j >= 0 and pattern[j] == sequence[i + j]:
            j -= 1
        if j < 0:
            patterns.append((i+1, i+len(pattern)))
        i += max(1, j - bad_char[ord(sequence[i + j])])
    if patterns:
        return patterns
    else:
        return -1

def KMP(pattern, sequence):
    range_lis = []
    c = 0
    M = len(pattern)
    N = len(sequence)

    # create lps[] that will hold the longest prefix suffix
    # values for pattern
    lps = [0]*M
    j = 0  # index for pat[]

    # Preprocess the pattern (calculate lps[] array) --> HOW??
    computeLPSArray(pattern, M, lps)
    i = 0  # index for txt[]
    while (N - i) >= (M - j):
        if pattern[j] == sequence[i]:
            i += 1
            j += 1

        if j == M:
            range_lis.append((i-j+1,i-j+M))
            c+=1
            j = lps[j-1]

        # mismatch after j matches
        elif i < N and pattern[j] != sequence[i]:
            # Do not match lps[0..lps[j-1]] characters,
            # they will match anyway
            if j != 0:
                j = lps[j-1]
            else:
                i += 1
    return range_lis

def computeLPSArray(pat, M, lps):
    len = 0  # length of the previous longest prefix suffix

    lps[0] = 0 # lps[0] is always 0
    i = 1

    # the loop calculates lps[i] for i = 1 to M-1
    while i < M:
        if pat[i] == pat[len]:
            len += 1
            lps[i] = len
            i += 1
        else:
            # This is tricky. Consider the example.
            # AAACAAAA and i = 7. The idea is similar
            # to search step.
            if len != 0:
                len = lps[len-1]

                # Also, note that we do not increment i here
            else:
                lps[i] = 0
                i += 1