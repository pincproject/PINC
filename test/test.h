/**
 * @file	    test.h
 * @author	    Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright   University of Oslo, Norway
 * @brief	    PINC unit testing framework
 * @date        16.12.15
 */

#ifndef TEST_H
#define TEST_H

#include <stdarg.h>

/**
 * @brief Tests an assertion.
 * @param	test		Boolean statement which should be true
 * @param	format		Printf-like format specifier for error message
 * @param	...			Printf-like arguments
 *
 * Unit tests are implemented as one function per test case. The test functions
 * should take no arguments and return 0 (int) if successful.
 *
 * Several assertions can be made per test case using this macro, which is
 * basically a wrapper for utAssertInner(). The developer makes an assertion by
 * calling this function with a test statement which should be true in case
 * everything works. If it _is_ true this macro does nothing. If it is false,
 * however, an error specified by the developer will be shown along with
 * debugging information. The error is specified in a printf-like way by using
 * a format specifier and arguments. An example test case is shown below.
 *
 * @code
 *	int myProduct(int a, int b){
 *		return a*b;
 *	}
 *
 *	static int testMyProduct(){
 *
 *		double result = myProduct(2,3);
 *		utAssert(result==6,"wrong result, 2*3 claimed to be %i.",result);
 *
 *		result = myProduct(30000,30000);
 *		utAssert(result==900000000,"doesn't work on large integers.");
 *
 *		return 0;
 *	}
 * @endcode
 *
 * The test is then performed by:
 *
 * @code
 *	utRun(&testMyProduct);
 * @endcode
 *
 * The macro will make the test function stop and return 1 on the first error.
 * Be careful when comparing floating point numbers as they may be slightly
 * different due to machine precision (machine epsilon) and finite accuracy in
 * numerical schemes. Rather than using a statement like "result==expected", use
 * "result-expected<precision".
 */
#define utAssert(...) do { int fail = utAssertInner(__FILE__,__func__,__LINE__,__VA_ARGS__); if(fail) return 1; } while (0)

/**
 * @brief Tests an assertion. Call indirectly using utAssert() macro.
 * @param	file	Name of file where assertion is (passed by macro)
 * @param	func	Name of function where assertion is (passed by macro)
 * @param	line	Line number where assertion is (passed by macro)
 * @param	test	Boolean statement which should be true
 * @param	format	Printf-like format specifier for error message
 * @param	...		Printf-like arguments
 * @return			0 if successful, 1 if failure
 * @see utAssert()
 */
int utAssertInner(const char *file, const char *func, int line, int test, const char * restrict format,...);

/**
 * @brief	Perform a test
 * @param	fun		Pointer to test function
 * @return			void
 *
 * A test is performed by calling this function with the address of the test
 * function as argument, rather than calling the test function directly.
 * This is used to make a summary of all tests.
 *
 * @see utAssert()
 */
void utRun(int (*fun)());

/**
 * @brief	Prints a summary of the tests
 * @return	void
 */
void utSummary();

/**
 * @brief	Performs all tests in aux.test.c
 * @return	void
 *
 * This prevents many small global test functions.
 */
void testAux();

/**
 * @brief	Performs all tests in grid.test.c
 * @return	void
 *
 * This prevents many small global test functions.
 */
void testGrid();

/**
 * @brief	Performs all tests in population.test.c
 * @return	void
 *
 * This prevents many small global test functions.
 */
void testPopulation();

/**
 * @brief	Performs all tests in io.test.c
 * @return	void
 *
 * This prevents many small global test functions.
 */
void testIo();

/**
 * @brief	Performs all tests in pusher.test.c
 * @return	void
 *
 * This prevents many small global test functions.
 */
void testPusher();

#endif // TEST_H
