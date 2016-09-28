
# HW1 Solution by BlackRice

## Group Info

| Group Member | SFU mail   | 
| ------ | ------ | 
| April Wang | ayw7@sfu.ca | 
| Jiahao Ke | jiahaok@sfu.ca |  
| Yu Tang | yta47@sfu.ca  |  
| Ruoxin Zhou | ruoxinz@sfu.ca  | 

## Solution

  In this project, we present a solution to segment Chinese words and our solution finally hits 89.747 on Leaderboard.
  
  We use some optimizations listed below to improve the naive segmentor based on the pseudo-code provided with HW 1 that uses unigram probabilities:
    
  - Using the bigram model to score word segmentation candidates.
  - Adding Good-Turing Algorithm to smooth probability data used by bigram and unigram model
  - Using Jelinek-Mercer Smoothing to improve the bigram model
  - Merging all the digit numbers appeared continuously in a line
  - Trying to recognize proper noun that is not collected in the dictionary

  For Jelinek-Mercer Smoothing, we write a bash script to test how different value of lambda would effect accuracy of final result. We get lambda = 0.17 for the best result.

  For merging all the digit numbers, we simply go through the lines that has been segmented and whenever digit numbers appearing continuously, we combine them together.

  For recognizing proper noun, we notice that the testing data is made of by pieces of news. For each news, there would likely be proper noun that appeared more than twice or three times. So for each segmented paragraph, we find those sets of words, while the length of each word is usually one but those sets of words appear continuously for more than twice. We consider those sets of words are likely to be proper noun and we combine them together.