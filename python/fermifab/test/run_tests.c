#include <stdbool.h>
#include <stdio.h>


typedef int (*test_function_t)(void);

// test function declarations
int TestTensorOp(void);


int main()
{
	test_function_t tests[] = {
		TestTensorOp
	};

	bool pass = true;

	size_t i;
	for (i = 0; i < sizeof(tests)/sizeof(tests[0]); i++)
	{
		int status = tests[i]();

		// print status message
		if (status < 0)
		{
			printf("Test crashed with internal problem!\n");
		}
		else if (status > 0)
		{
			printf("Test failed!\n");
		}
		else    // status == 0
		{
			printf("Test passed.\n");
		}

		printf("________________________________________________________________________________\n");

		if (status != 0) {
			pass = false;
		}
	}

	if (pass)
	{
		printf("All tests passed :-)\n");
	}
	else
	{
		printf("At least one test failed or crashed :-(\n");
	}

	return 0;
}
