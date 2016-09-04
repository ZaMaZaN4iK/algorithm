/*
  Copyright (c) Alexander Zaitsev <zamazan4ik@gmail.com>, 2016
  Distributed under the Boost Software License, Version 1.0. (See
  accompanying file LICENSE_1_0.txt or copy at
  http://www.boost.org/LICENSE_1_0.txt)
  See http://www.boost.org/ for latest version.
*/

#ifndef BOOST_ALGORITHM_SUFFIX_AUTOMATA_HPP
#define BOOST_ALGORITHM_SUFFIX_AUTOMATA_HPP

#include <map>
#include <unordered_map>
#include <memory>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/make_unique.hpp>

namespace boost { namespace algorithm {


template <typename T, template<typename ...> class Container, typename ...Args>
class suffix_automata
{
private:
    struct node
    {
        size_t length_;
        node* link_;
        Container<T, node*, Args...> next;

        node(size_t length = 0, node* link = nullptr) : length_(length), link_(link) {}
    };

public:
    using value_type = T;
    using node_type  = node;

    suffix_automata() : root(new node_type)
    {
    }

    template<typename ForwardIterator>
    explicit suffix_automata(ForwardIterator begin, ForwardIterator end) : suffix_automata()
    {
        init(begin, end);
    }

    template<typename Range>
    explicit suffix_automata(const Range& range) :
            suffix_automata(boost::begin(range), boost::end(range))
    {}

    template<typename ForwardIterator>
    void init(ForwardIterator begin, ForwardIterator end)
    {
        last = root;
        while(begin != end)
        {
            //cur = sz++;
            cur = new node_type;
            cur->length_ = last->length_ + 1;
            value_type& c = *begin;
            for(p = last; p && !p->next.count(c); p = p->link_)
            {
                p->next[c] = cur;
            }
            if(!p)
            {
                cur->link_ = root;
            }
            else
            {
                node_type* q = p->next[c];
                if(p->length_ + 1 == q->length_)
                {
                    cur->link_ = q;
                }
                else
                {
                    //ll clone = sz++;
                    node_type* clone = new node_type;
                    clone->length_  = p->length_ + 1;
                    clone->link_ = q->link_;
                    clone->next = q->next;
                    for(;p && p->next[c] == q; p = p->link_)
                    {
                        p->next[c] = clone;
                    }
                    cur->link_ = q->link_ = clone;
                }
            }
            last = cur;
            ++begin;
        }
    }

    template<typename ForwardIterator>
    bool find(ForwardIterator begin, ForwardIterator end)
    {
        node_type* v = root;
        while(begin != end)
        {
            if(!v->next.count(*begin))
            {
                return false;
            }
            else
            {
                v = v->next[*begin];
            }
            ++begin;
        }
        return true;
    }

    template<typename Range>
    bool find(const Range& range)
    {
        return find(boost::begin(range), boost::end(range));
    }

private:
    node_type *cur, *last, *p, *root;
};


//Object interface
template <typename T, typename Pred = std::less<T>>
using suffix_automata_map = suffix_automata<T, std::map, Pred>;

template <typename T, typename Hash = std::hash<T>, typename Comp = std::equal_to<T>>
using suffix_automata_hashmap = suffix_automata<T, std::unordered_map, Hash, Comp>;


//Functional interface

/// \fn aho_corasick_map ( RAIterator corpus_begin, RAIterator corpus_end,
///                        ForwardIterator pat_begin, ForwardIterator pat_end,
///                        Callback cb)
/// \return true if all callback callings return true, else false
///
/// \param corpus_begin The start of the corpus sequence
/// \param corpus_end   One past the end of the corpus sequence
/// \param pat_begin	The start of the patterns sequence
/// \param pat_end	    One past the end of the patterns sequence
/// \param cb 		    Callback for matches
///
template <typename T, typename Predicate = std::less<T>, typename RAIterator,
        typename ForwardIterator, typename Callback>
bool suffix_automata_map ( RAIterator corpus_begin, RAIterator corpus_end,
                        ForwardIterator pat_begin, ForwardIterator pat_end,
                        Callback cb)
{
    suffix_automata<T, std::map, Predicate> obj(pat_begin, pat_end);
    return obj.find(corpus_begin, corpus_end, cb);
}

/// \fn aho_corasick_hashmap ( RAIterator corpus_begin, RAIterator corpus_end,
///                            ForwardIterator pat_begin, ForwardIterator pat_end,
///                            Callback cb)
/// \return true if all callback callings return true, else false
///
/// \param corpus_begin The start of the corpus sequence
/// \param corpus_end   One past the end of the corpus sequence
/// \param pat_begin	The start of the patterns sequence
/// \param pat_end	    One past the end of the patterns sequence
/// \param cb 		    Callback for matches
///
template <typename T, typename Hash = std::hash<T>, typename Comp = std::equal_to<T>, typename RAIterator,
        typename ForwardIterator, typename Callback>
bool aho_corasick_hashmap ( RAIterator corpus_begin, RAIterator corpus_end,
                            ForwardIterator pat_begin, ForwardIterator pat_end,
                            Callback cb)
{
    AhoCorasick<T, std::unordered_map, Hash, Comp> obj(pat_begin, pat_end);
    return obj.find(corpus_begin, corpus_end, cb);
}

}}

#endif //BOOST_ALGORITHM_SUFFIX_AUTOMATA_HPP
