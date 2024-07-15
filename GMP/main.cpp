#include <iostream>
#include <fstream>
#include <future>
#include <gmp.h>
#include <chrono>
#include <math.h>
#include <string.h>
// 1000000 digits -> 75864 iterations
const unsigned long DIGITS = 1000000; // 250000 took 1938435ms vs 1180054ms
const unsigned long ITERATIONS = (unsigned long)ceil(DIGITS / log10(151931373056000));
const unsigned long PERCISION = DIGITS*log2(10);

// https://stackoverflow.com/questions/7459259/inserting-characters-into-a-string
void append(char subject[], const char insert[], int pos) {
    char *buf = (char*)malloc(strlen(subject) + strlen(insert) + 2);
    memset(buf, 0, 100);

    strncpy(buf, subject, pos); // copy at most first pos characters
    int len = strlen(buf);
    strcpy(buf+len, insert); // copy all of insert[] at the end
    len += strlen(insert);  // increase the length by length of insert[]
    strcpy(buf+len, subject+pos); // copy the rest

    strcpy(subject, buf);   // copy it back to subject
    // Note that subject[] must be big enough, or else segfault.
    // deallocate buf[] here, if used malloc()
    // e.g. free(buf);
}

void compute_sqrt_part(mpf_t* result)
{
	mpf_t a,b;
	mpf_init(b);
	mpf_init_set_ui(a, 426880);

	// sqrt(10005)
	mpf_sqrt_ui(b, 10005);

	// 426880 * sqrt(10005)
	mpf_mul(*result, a, b);

	gmp_printf("DONE: Computed sqrt part.\n");
}

void compute_P(mpz_t a, mpz_t b, mpz_t* result)
{
	// temp
	mpz_t c, d, e;
	mpz_init_set_ui(c, 0);
	mpz_init_set_ui(d, 0);
	mpz_init_set_ui(e, 0);

	// set result to 1 or the result is always 0
	mpz_set_ui(*result, 1);

	// j = a
	mpz_t j;
	for(mpz_init_set(j, a);mpz_cmp(j, b)<0;mpz_add_ui(j, j, 1))
	{
		// -(6j-1)
		mpz_set_si(d, -1);
		mpz_mul_ui(c, j, 6);
		mpz_sub_ui(c, c, 1);
		
		mpz_mul(c, c, d);
		// (2j-1)
		mpz_mul_ui(d, j, 2);
		mpz_sub_ui(d, d, 1);
		// (6j-5)
		mpz_mul_ui(e, j, 6);
		mpz_sub_ui(e, e, 5);

		// (6j-1)(2j-1)(6j-5)
		mpz_mul(c, c, d);
		mpz_mul(c, c, e);

		// (6j-1)(2j-1)(6j-5)

		// result *= -(6j-1)(2j-1)(6j-5)
		mpz_mul(*result, *result, c);
	}

	mpz_clear(c);
	mpz_clear(d);
	mpz_clear(e);
	mpz_clear(j);
}

void compute_Q(mpz_t a, mpz_t b, mpz_t* result)
{
	// temp
	mpz_t c;
	mpz_init_set_ui(c, 0);

	// set result to 1 or the result is always 0
	mpz_set_ui(*result, 1);

	// https://www.wolframalpha.com/input?i=divisors+of+10939058860032000
	// set the 10939058860032000 to the split version
	mpz_t _10939058860032000;
	mpz_init_set_ui(_10939058860032000, 32768);
	mpz_mul_ui(_10939058860032000, _10939058860032000, 9);
	mpz_mul_ui(_10939058860032000, _10939058860032000, 125);
	mpz_mul_ui(_10939058860032000, _10939058860032000, 12167);
	mpz_mul_ui(_10939058860032000, _10939058860032000, 24389);

	// j = a
	mpz_t j;

	for(mpz_init_set(j, a);mpz_cmp(j, b)<0;mpz_add_ui(j, j, 1))
	{
		// j ^ 3
		mpz_mul(c, j, j);
		mpz_mul(c, c, j);

		mpz_mul(c, c, _10939058860032000);
		
		// result *= c
		mpz_mul(*result, *result, c);
	}

	mpz_clear(c);
	mpz_clear(_10939058860032000);
	mpz_clear(j);
}

void compute_S(mpz_t a, mpz_t b, mpf_t* result)
{
	// temp
	mpf_t c, d, e;
	mpf_init_set_ui(c, 0);
	mpf_init_set_ui(d, 0);
	mpf_init_set_ui(e, 0);

	mpz_t f;
	mpz_init_set_ui(f, 0);

	mpz_t p, q;
	mpz_init(p);
	mpz_init(q);

	// j = a
	mpz_t j;

	for(mpz_init_set(j, a);mpz_cmp(j, b)<0;mpz_add_ui(j, j, 1))
	{
		// j + 1
		mpz_set(f, j);
		mpz_add_ui(f, f, 1);

		std::future<void> p_task = std::async(std::launch::async, compute_P, a, f, &p);
		std::future<void> q_task = std::async(std::launch::async, compute_Q, a, f, &q);

		p_task.get();
		q_task.get();
		
		// 545140134k
		mpz_mul_ui(f, j, 545140134);
		// 545140134k + 13591409
		mpz_add_ui(f, f, 13591409);

		// cast to float
		mpf_set_z(c, p);
		mpf_set_z(d, q);
		mpf_set_z(e, f);

		// P / Q
		mpf_div(c, c, d);
		// (P / Q) * (545140134k + 13591409)
		mpf_mul(c, c, e);

		// result += (P / Q) * (545140134k + 13591409)
		mpf_add(*result, *result, c);

		gmp_printf("S iteration: %Zd\n", j);
	}

	mpf_clear(c);
	mpf_clear(d);
	mpf_clear(e);

	mpz_clear(j);
	mpz_clear(f);
	mpz_clear(p);
	mpz_clear(q);
}

void compute(mpf_t& pi, mpz_t n)
{
	mpz_t one;
	mpz_init_set_ui(one, 1);

	mpf_t sqrt_part;
	mpf_init(sqrt_part);

	std::future<void> sqrt_part_task = std::async(std::launch::async, compute_sqrt_part, &sqrt_part);

	mpf_t s;
	mpf_init(s);
	
	std::future<void> s_task = std::async(std::launch::async, compute_S, one, n, &s);
	s_task.get();

	// (13591409 + S)
	mpf_add_ui(s, s, 13591409);
	
	sqrt_part_task.get();

	// pi = (426880 * sqrt(10005)) / (13591409 + S)
	mpf_div(pi, sqrt_part, s);
}

int main ()
{
	mpf_set_default_prec(PERCISION);

	std::cout << mpf_get_default_prec() << std::endl;

	mpf_t result;
	mpf_init_set_ui(result, 0);

	mpz_t iterations;
	mpz_init_set_ui(iterations, ITERATIONS);

	gmp_printf("Running %Zd iterations. \n", iterations);

	auto start = std::chrono::high_resolution_clock::now();

	compute(result, iterations);

	auto end = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	auto millis = duration.count();
	auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
	auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
	auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
	gmp_printf("Computation of %Zd iterations took %dms(%ds)(%dmin)(%dh).", iterations, millis, seconds, minutes, hours);

	mp_exp_t exponent;
	// Convert to string
	char* pi_string = mpf_get_str(NULL, &exponent, 10, 0, result);
	// Insert . to the second place
	append(pi_string, ".", 1);

	std::ofstream output_file;
	output_file.open("pi.txt");


	output_file << pi_string;

	output_file.close();

	mpf_clear(result);
	free(pi_string);

	return 0;
}
