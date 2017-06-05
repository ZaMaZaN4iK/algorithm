/*
   Copyright (c) Evgenii Zheltonozhskii 2017.

   Distributed under the Boost Software License, Version 1.0. (See accompanying
   file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

    For more information, see http://www.boost.org
*/

#ifndef BOOST_ALGORITHM_FUZZY_SEARCH_H
#define BOOST_ALGORITHM_FUZZY_SEARCH_H

#include <algorithm>
#include <functional>
#include <map>
#include <unordered_set>
#include <vector>

/**
 *
 * @tparam ForwardIt1
 * @tparam SearchedType
 * @param first
 * @param last
 * @param k
 * @param q
 * @param i
 * @param s
 * @param S1
 * @param M
 * @param Ds
 * @param D
 */
template <class ForwardIt1, class SearchedType>
void preprocess_helper_mismatch(ForwardIt1 first,
                                ForwardIt1 last,
                                size_t m,
                                size_t k,
                                size_t q,
                                size_t i,
                                std::vector<SearchedType> &s,
                                std::unordered_set<SearchedType> &S1,
                                std::map<std::vector<SearchedType>, size_t> &M,
                                std::map<std::vector<SearchedType>, size_t> &Ds,
                                std::vector<std::vector<size_t>> &D) {
  if (i == q + 1) {
    M[s]  = D[q][m];
    Ds[s] = std::find_if(D[q].rbegin() + 1, D[q].rend(), [k](size_t j) { return j <= k; }) - D[q].rbegin();
  } else {
    s.push_back(*S1.begin());
    for (auto &c : S1) {
      for (int j = 1; j <= m; ++j) {
        D[i][j] = D[i - 1][j - 1] + ((c == *std::next(first, j - 1)) ? 0 : 1); //TODO: 1D
      }
      // TODO: For maximum performance it is crucial how the value of a q-gram is computed during searching.
      s.back() = c;
      preprocess_helper_mismatch(first, last, m, k, q, i + 1, s, S1, M, Ds, D);
    }
    s.pop_back();
  }
}

/**
 *
 * @tparam ForwardIt1
 * @tparam SearchedType
 * @param first
 * @param last
 * @param m
 * @param k
 * @param q
 * @param i
 * @param s
 * @param S1
 * @param M
 * @param Ds
 * @param D
 */
template <class ForwardIt1, class SearchedType>
void preprocess_helper_difference(ForwardIt1 first,
                                  ForwardIt1 last,
                                  size_t m,
                                  size_t k,
                                  size_t q,
                                  size_t i,
                                  std::vector<SearchedType> &s,
                                  std::unordered_set<SearchedType> &S1,
                                  std::map<std::vector<SearchedType>, size_t> &M,
                                  std::map<std::vector<SearchedType>, size_t> &Ds,
                                  std::vector<std::vector<size_t>> &D) {
  if (i == q + 1) {
    M[s]  = D[q][m];
    Ds[s] = std::find_if(D[q].rbegin() + 1, D[q].rend(), [k](size_t j) { return j <= k; }) - D[q].rbegin();
  } else {
    s.push_back(*S1.begin());
    for (auto &c : S1) {
      for (int j = 1; j <= m; ++j) {
        size_t a1 = D[i - 1][j] + 1, a2 = D[i][j - 1] + 1,
            a3 = D[i - 1][j - 1] + ((c == *std::next(first, j - 1)) ? 0 : 1);
        D[i][j]   = std::min({a1, a2, a3});
      }
      // TODO: For maximum performance it is crucial how the value of a q-gram is computed during searching.
      s.back() = c;
      preprocess_helper_difference(first, last, m, k, q, i + 1, s, S1, M, Ds, D);
    }
    s.pop_back();
  }
}
/**
 *
 * @tparam ForwardIt1
 * @tparam SearchedType
 * @param first
 * @param last
 * @param k
 * @param q
 * @param S1
 * @param M
 * @param Ds
 */
template <class ForwardIt1, class SearchedType>
void preprocess_mismatch(ForwardIt1 first,
                         ForwardIt1 last,
                         size_t k,
                         size_t q,
                         std::unordered_set<SearchedType> &S1,
                         std::map<std::vector<SearchedType>, size_t> &M,
                         std::map<std::vector<SearchedType>, size_t> &Ds) {
  size_t m = std::distance(first, last);
  std::vector<std::vector<size_t>> D(q + 1, std::vector<size_t>(m + 1, 0));
  auto s = std::vector<SearchedType>();
  s.reserve(q);
  preprocess_helper_mismatch(first, last, m, k, q, 1, s, S1, M, Ds, D);
}

/**
 *
 * @tparam ForwardIt1
 * @tparam SearchedType
 * @param first
 * @param last
 * @param k
 * @param q
 * @param S1
 * @param M
 * @param Ds
 */
template <class ForwardIt1, class SearchedType>
void preprocess_difference(ForwardIt1 first,
                           ForwardIt1 last,
                           size_t k,
                           size_t q,
                           std::unordered_set<SearchedType> &S1,
                           std::map<std::vector<SearchedType>, size_t> &M,
                           std::map<std::vector<SearchedType>, size_t> &Ds) {
  size_t m = std::distance(first, last);
  std::vector<std::vector<size_t>> D(q + 1, std::vector<size_t>(m + 1, 0));

  auto s = std::vector<SearchedType>();
  s.reserve(q);
  preprocess_helper_difference(first, last, m, k, q, 1, s, S1, M, Ds, D);
}

/**
 *
 * @tparam ForwardIt1
 * @tparam ForwardIt2
 * @param first
 * @param s_first
 * @param k
 * @param m
 * @return
 */
template <class ForwardIt1, class ForwardIt2>
bool validate_mismatch(ForwardIt1 first, ForwardIt2 s_first, size_t k, size_t m) {  // TODO bidir specialization
  auto c = 0;
  for (int i = m - 1; i >= 0; --i) {
    if (*std::next(first, i) != *std::next(s_first, i)) {  // TODO
      ++c;
      if (c > k) return false;
    }
  }
  return true;
}

/**
 *
 * @tparam ForwardIt1
 * @tparam ForwardIt2
 * @param first
 * @param s_first
 * @param k
 * @param m
 * @return
 */
template <class ForwardIt1, class ForwardIt2>
bool validate_difference(ForwardIt1 first, ForwardIt2 s_first, size_t k, size_t m) {  // TODO bidir specialization
  std::vector<size_t> line1(m + 1, 0), line2(m + 1, 0);
  std::vector<size_t> *curr=&line1, *prev=&line2;
  for (int i = 1; i < m + k + 1; ++i) {
    for (int j = 1; j < m + 1; ++j) {
      size_t a1 = (*prev)[j] + 1, a2 = (*curr)[j - 1] + 1,
          a3 = (*prev)[j - 1] + ((*std::next(first, i - 1) == *std::next(s_first, j - 1)) ? 0 : 1);
      (*curr)[j]   = std::min({a1, a2, a3});
    }
    std::swap(curr, prev);
  }
  return (*prev)[m] <= k;  // TODO
}

/**
 *
 * @tparam ForwardIt1
 * @tparam ForwardIt2
 * @param first The start of the data to search in
 * @param last The end of the data to search in
 * @param s_first The start of the data to search
 * @param s_last The end of the data to search
 * @param k possible number of mistakes
 * @return
 */
template <class ForwardIt1, class ForwardIt2>
ForwardIt1 fuzzy_search(
    ForwardIt1 first, ForwardIt1 last, ForwardIt2 s_first, ForwardIt2 s_last, size_t k, bool mismatch) {
  size_t m = std::distance(s_first, s_last);  // length of needle

  using SearchedType  = typename std::iterator_traits<ForwardIt2>::value_type;
  using SearchingType = typename std::iterator_traits<ForwardIt1>::value_type;
  auto preprocess =
      mismatch ? preprocess_mismatch<ForwardIt2, SearchedType> : preprocess_difference<ForwardIt2, SearchedType>;
  auto validate = mismatch ? validate_mismatch<ForwardIt1, ForwardIt2> : validate_difference<ForwardIt1, ForwardIt2>;

  // TODO optimal parameters and corner cases
  size_t reduced_alphabet_size = 16;  // size of reduced alphabet
  size_t q                     = 5;   // size of qgrams
  // should depend on searched data size

  size_t start_search_position = mismatch ? m - q : m - k - q;

  std::unordered_set<SearchedType> S(s_first, s_last);  // S is set of all chars in P
  std::vector<SearchingType> T1(first, last);           // T1 is copy of T

  std::map<SearchedType, size_t> count;
  bool added = false;
  SearchingType extra;  // extra character replacing everything not present in P
  for (auto it = T1.begin(); it != T1.end(); ++it) {
    if (S.find(*it) == S.end()) {
      if (added) {
        *it = extra;
      } else {
        added = true;
        extra = *it;
        S.insert(*it);  // add extra character to S. Now S is pattern alphabet
      }
    }
    count[*it] += 1;
  }

  std::vector<std::pair<size_t, size_t>> freq(count.begin(), count.end());
  // sort by frequency
  std::sort(freq.begin(), freq.end(),
            [](std::pair<size_t, size_t> a, std::pair<size_t, size_t> b) { return b.second < a.second; });

  std::unordered_set<SearchedType> S1;           // S1 is reduced pattern alphabet (Sigma')
  std::map<SearchedType, SearchedType> mapping;  //  mapping  Sigma -> Sigma'
  for (size_t i = 0; i < reduced_alphabet_size; ++i) {
    mapping[freq[i].first] = freq[i].first;  // most frequent chars are mapped to themselves
    S1.insert(freq[i].first);
  }
  std::vector<std::pair<size_t, size_t>> freq_reduced(
      freq.begin(),
      freq.begin() + reduced_alphabet_size);  // frequency of reduced alphabet
  for (size_t i = reduced_alphabet_size; i < freq.size(); ++i) {
    mapping[freq[i].first] = freq_reduced.back().first;  // most frequent char mapped to least frequent in reduced
    freq_reduced.back().second += freq[i].second;
    std::sort(freq_reduced.begin(), freq_reduced.end(), [](std::pair<size_t, size_t> a, std::pair<size_t, size_t> b) {
      return b.second < a.second;
    });  // TODO: optimal?
  }

  std::vector<SearchedType> P1(s_first, s_last);  // P1 is copy of P1
  for (auto it = P1.begin(); it != P1.end(); ++it) {
    *it = mapping[*it];
  }

  std::map<std::vector<SearchedType>, size_t> M, Ds;
  preprocess(s_first, s_last, k, q, S1, M, Ds);

  for (auto it = T1.begin(); it != T1.end(); ++it) {
    *it = mapping[*it];
  }

  auto s = std::next(T1.begin(), start_search_position);
  while (std::distance(s, T1.end()) > q) {
    std::vector<SearchingType> current_gram(s, std::next(s, q));
    if (M[current_gram] <= k) {
      auto orig_s = std::next(first, std::distance(T1.begin(), s + q - m - k));
      if (validate(orig_s, s_first, k, m)) return orig_s;
    }
    if (std::distance(s, T1.end()) <= Ds[current_gram] + q - 1) break;
    s = std::next(s, Ds[current_gram]);
  }

  return last;
}

#endif  // BOOST_ALGORITHM_FUZZY_SEARCH_H