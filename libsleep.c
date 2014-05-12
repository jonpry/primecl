#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* 
 * To compile run:
 * gcc -O2 -fPIC -shared -Wl,-soname,libsleep.so -o libsleep.so libsleep.c
 * 
 * To use:
 * LD_PRELOAD="./libsleep.so" ./cgminer
 * 
 * You can configure sleep time by setting 
 * YIELD_SLEEP_TIME environment variable (in microseconds)
 * Default is 1000usec
 * Example:
 * YIELD_SLEEP_TIME="1500" LD_PRELOAD="./libsleep.so" ./cgminer
 *
 * Tips are welcome: 
 * Bitcoin:  1FQMFpqnCH1ATPGoHLnmSAXB1gBh9vAKXC 
 * Litecoin: LZiXcRvUr5wJrkfvc9vXvqDwe51YfVPRkv
 *
 */

useconds_t yield_sleep_time = 1000;

static void __attribute__ ((constructor)) lib_init(void) {
	int stime = 0;
	char *env_stime = getenv("YIELD_SLEEP_TIME");

	if(env_stime) {
		stime = atoi(env_stime);

		if(stime > 0)
			yield_sleep_time = stime;
	}

	printf("libsleep: Sleep time: %uusec\n", yield_sleep_time);
}

int sched_yield(void) {
	usleep(yield_sleep_time);
	return 0;
}
