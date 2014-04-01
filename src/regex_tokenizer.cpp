#include "regex_tokenizer.hpp"

void regex_tokenizer::set_regex( std::string exp, bool case_sensitive )
{
        try
        {
                if(!case_sensitive)
                        re.assign(exp,boost::regex_constants::icase);
                else
                        re.assign(exp);
        }
        catch ( ... )
        {
                BOOST_THROW_EXCEPTION(module_not_found() 
                                  << errstr_info( std::string("Invalid regular expression ") + exp)
                                    );
        }

}

std::string regex_tokenizer::get_regex()
{
    return re.expression();
}



regex_tokenizer::regex_tokenizer() :
        FLOAT_REGEX("-?[0-9]+\\.?[0-9]*")
{

}

regex_tokenizer::regex_tokenizer( std::string regex, bool case_senstive /*= true*/ ) :
        FLOAT_REGEX("-?[0-9]+\\.?[0-9]*")
{
        try
        {
                set_regex(regex,case_senstive);
        }
        catch ( ... )
        {
                //nothing. don't want a zombie
        }

}
regex_tokenizer::~regex_tokenizer()
{

}
