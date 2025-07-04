#include <gtest/gtest.h>
#include "daily_accumulator.hpp"

class data
{
public:
	bool is_new_day();
	int steps_per_day();

	void set_steps(int& step)
	{ _steps_per_day = step; };
	
	int get_steps()
	{ return _steps_per_day; };

	void step_forward() 
	{ ++_count; };
	
	void reset()
	{ _count = 0; };

private:
	int _steps_per_day;
	int _count = 0;
};

int data::steps_per_day()
{
    return get_steps();
};

bool data::is_new_day()
{
	return _count % _steps_per_day == 0 && _count > 0;
};

class DailyAccumulatorTest : public testing::Test
{
protected:
    template<typename T>
	double sum_to_step(int i,int steps_per_day,T arr)
	{
		int start = nearest_day_start(i,d.get_steps());
		double sum = 0;
		for (int j = start; j < i; ++j)
		{
			sum += arr.at(j);
		}
		return sum;
	};

	int nearest_day_start(int i, int step_per_day)
	{
		if (step_per_day == 0) {return 0;}

		int day_start_index = (i / step_per_day) * step_per_day;

		return (day_start_index >= step_per_day) ? day_start_index : 0;
	};

    data d;
};

TEST_F(DailyAccumulatorTest,HelperFunctionSumToStepTest)
{
    EXPECT_EQ(1,1);
	const int steps = 5;

	std::vector<double> T{-1.0,0.0,2.0,3.0,6.2};

	// sum to step test
	double sum = sum_to_step(0,steps,T);

	ASSERT_EQ(sum,0.0);

	sum = sum_to_step(2,steps,T);

	ASSERT_EQ(sum,-1.0);

	sum = sum_to_step(5,steps,T);

	ASSERT_EQ(sum,10.2);
};

TEST_F(DailyAccumulatorTest,HelperFunctionNearestDayStartTest)
{
	int start = nearest_day_start(3,12);

	ASSERT_EQ(start,0);

	start = nearest_day_start(99,400);

	ASSERT_EQ(start,0);

	start = nearest_day_start(29,12);

	ASSERT_EQ(start,24);

	start = nearest_day_start(35,4);

	ASSERT_EQ(start,32);

	start = nearest_day_start(35+4,4);

	ASSERT_EQ(start,36);
};

TEST_F(DailyAccumulatorTest,Test)
{
	double temperature = 0.0;
	data d;
	const int steps = 4;
    int steps_copy = steps;
	d.set_steps(steps_copy);

	daily_accumulator<data> accumulate_temperature(d);
	accumulate_temperature.bind_to_var(temperature);

	std::array<double,2*steps> T{-1.0,-3.2,4.0,1.2,-2.9,-10.1,-5.1,5.0};

	for (int i = 0; i < 2 * steps; ++i)
	{
		temperature = T[i];
		accumulate_temperature.execute();
		d.step_forward();
        std::cout << "test: " << i % steps << "," << i << std::endl;
		if (i < steps)
        {
            ASSERT_EQ(accumulate_temperature.get_last_mean(),0.0) << "Loop: " << i;
        }
        else if (i % steps == 0)
		{
            std::cout << "har" <<std::endl;
            //on new day, see if sum computed
			double sum_total = sum_to_step(i,steps,T);
			ASSERT_DOUBLE_EQ(sum_total / steps,accumulate_temperature.get_last_mean()) << "Loop: " << i;
		}
		//else
		//{
        //    std::cout << "here" << std::endl;
	//		// Length of array fixed at 2*steps
	//		double sum_total = (T[0] + T[1] + T[2] + T[3])/steps;
//			ASSERT_EQ(accumulate_temperature.get_last_mean(),sum_total) << "Loop: " << i;
//		}
	}
};
