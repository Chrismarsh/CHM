#pragma once

#include <string>
#include <vector>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/any.hpp>

#include "exception.hpp"

// Class: regex_tokenizer
//	A class that allows one to tokenize a string based on a regular expression
class regex_tokenizer
{
        public:
                const std::string FLOAT_REGEX;
                /*
                        Function: regex_tokenizer
                                 Creates an empty tokenizer

                        Parameters: 
                                None

                        Throws:
                                Never

                        Returns:   
                                 - 
                */
                regex_tokenizer();

                ~regex_tokenizer();


                /*
                        Function: regex_tokenizer
                                 Creates a tokenizer with a given regex

                        Parameters: 
                                std::string regex - A valid regular expression
                                bool case_senstive - default = true

                        Throws:
                                If an incorrect regex is entered, it will not throw. When you use <tokenize> an exception will be thrown.

                        Returns:   
                                 - 
                */
                regex_tokenizer(std::string regex, bool case_senstive = true);

                /*
                        Function: set_regex
                                Sets the regular expression to use. Does not throw on an invalid regex.

                        Parameters:
                                exp - A valid regular expression.
                                case_sensitive - Defines if a regex search should be case sensitive. Defaults to true (case sensitive).

                        Throws:
                                std::runtime_error on invalid regex

                        Returns:
                                void
                 */
                void set_regex(std::string exp,bool case_sensitive = true);


                /*
                        Function: get_regex
                                Returns the current regular expression

                        Parameters:
                                None.

                        Throws:
                                Never.

                        Returns:
                                std::string - String representation of the regular expression
                */
                std::string get_regex();


                /*
                        Function: tokenize
                                Tokenizes the given item. Each match is put into the returning vector

                        Parameters:
                                std::string item - String to tokenize.

                        Throws:
                                std::runtime_error - Thrown on either:
                                                                                        a) an invalid cast b/w type T and the result from the regular
                                                                                                expression match ( "asdf" -> int).
                                                                                        b) an invalid regular expression

                        Returns:
                                 A vector holding each match.
                */
                template<typename T>
                std::vector<T> tokenize( std::string item )
                {
                        boost::match_results<std::string::const_iterator> what;

                        std::string::const_iterator start = item.begin();
                        std::string::const_iterator end = item.end();

                        std::vector<T> items;



                        while (boost::regex_search(start, end, what, re))
                        {
                                try
                                {
                                    std::string tmp = what[0];
                                    items.push_back(boost::lexical_cast<T>(tmp));
                                }
                                catch (boost::bad_lexical_cast)
                                {
                                    BOOST_THROW_EXCEPTION(bad_lexical_cast() 
                                        << errstr_info( std::string("Bad lexical cast"))
                                            );
                                }


                                // Update the beginning of the range to the character
                                // following the match
                                start = what[0].second;
                        }

                        return items;
                }


                //TODO:Copy constructor

private:
                //regular expression
                boost::regex re;

};
