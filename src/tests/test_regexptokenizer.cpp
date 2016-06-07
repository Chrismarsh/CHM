
#include "regex_tokenizer.hpp"
#include "logger.hpp"
#include "gtest/gtest.h"

class RegexpTokenizerTest : public testing::Test
{
protected:

    virtual void SetUp()
    {
        logging::core::get()->set_logging_enabled(false);
    }

};

TEST_F(RegexpTokenizerTest,InavildRegEx)
{
    std::string bad_text(1024, ' ');

    regex_tokenizer token;

    token.set_regex("(.+)+xyz");

    ASSERT_ANY_THROW(token.tokenize<std::string>(bad_text));
}

TEST_F(RegexpTokenizerTest,TestInt)
{
    regex_tokenizer token;

    std::string line1 = "1 2	3asd		4";
    std::vector<int> items;

    token.set_regex("[0-9]+");
    items = token.tokenize<int>(line1);

    ASSERT_EQ(items.size() , 4);
    ASSERT_EQ(items[0] , 1);
    ASSERT_EQ(items[1] , 2);
    ASSERT_EQ(items[2] , 3);
    ASSERT_EQ(items[3] , 4);
}

TEST_F(RegexpTokenizerTest,BadCast)
{
    regex_tokenizer token;

    std::string line1 = "filename.cpp asdf";
    std::vector<int> items;

    token.set_regex("[a-zA-Z]+");
    ASSERT_ANY_THROW(items = token.tokenize<int>(line1));

    ASSERT_EQ(items.size() , 0);
}
TEST_F(RegexpTokenizerTest,TestDouble)
{
    regex_tokenizer token;

    std::string line1 = "1.1 2.2	3.3333asd		4.4";
    std::vector<double> items;
    items.clear();

    token.set_regex("[0-9]+\\.[0-9]+");
    items = token.tokenize<double>(line1);

    ASSERT_EQ(items.size() , 4);
    ASSERT_DOUBLE_EQ(items[0] , 1.1);
    ASSERT_DOUBLE_EQ(items[1] , 2.2);
    ASSERT_DOUBLE_EQ(items[2] , 3.3333);
    ASSERT_DOUBLE_EQ(items[3] , 4.4);
}