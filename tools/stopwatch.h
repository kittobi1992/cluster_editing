#if defined(__unix) && defined(__sun)
#include <sys/time.h>
#elif defined(__linux__)
#include <sys/time.h>
#include <sys/resource.h>
#endif


/**
 * Time measurement for real ("wall-clock") time. As is the nature of a
 * stopwatch, only time @b differences can be measured, no absolute times.
 *
 * Notes (valid also for ProcessStopwatch):
 * - Implementation on Solaris is via gethrtime() and gethrvtime()
 * - on Linux: gettimeofday() and getrusage()
 * - gettimeofday() is not necessarily linear, it can even decrease under some
 *   circumstances. If that happens, elapsed() returns an incorrect (perhaps
 *   even negative) time.
 *
 * @author Marcel Martin <Marcel.Martin@CeBiTec.Uni-Bielefeld.DE>
 *
 * @ingroup utils
 */
class Stopwatch {
	public:
		/**
		 * Creates a new stopwatch. The stopwatch is automatically started,
		 * so there is no need to call start() directly after construction.
		 */
		Stopwatch() { start(); }
		void start();
		double elapsed() const;

	private:
#if defined(__unix) && defined(__sun)
		hrtime_t start_time;
#elif defined(__linux__)
		timeval tv_start;
#endif
};


/**
 * Starts the stopwatch.
 */
inline void Stopwatch::start() {
#if defined(__unix) && defined(__sun)
	start_time = gethrtime();
#elif defined(__linux__)
	gettimeofday(&tv_start, NULL);
#endif
}


/**
 * Measures real time elapsed between the last call to start() and the current time.
 * @return time elapsed in seconds.
 */
inline double Stopwatch::elapsed() const {
#if defined(__unix) && defined(__sun)
	hrtime_t now = gethrtime();
	// convert to seconds since gethrtime() returns nanoseconds
	return (now-start_time)/1e9;
#elif defined(__linux__)
	timeval tv_end;
	gettimeofday(&tv_end, NULL);
	timeval tv_elapsed;
	timersub(&tv_end, &tv_start, &tv_elapsed);
	return double(tv_elapsed.tv_sec) + tv_elapsed.tv_usec * 1e-6;
#endif
}


/**
 * Time measurements for process time (CPU time). Measures only time
 * @b differences, just like Stopwatch.
 *
 * @author Marcel Martin <Marcel.Martin@CeBiTec.Uni-Bielefeld.DE>
 * @ingroup utils
 */
class ProcessStopwatch {
	public:
		ProcessStopwatch() { start(); }
		void start();
		double elapsed() const;

	private:
#if defined(__unix) && defined(__sun)
		hrtime_t start_time;
#elif defined(__linux__)
		timeval tv_start;
#endif
};


/**
 * Starts the stopwatch.
 */
inline void ProcessStopwatch::start() {
#if defined(__unix) && defined(__sun)
	start_time = gethrvtime();
#elif defined(__linux__)
	rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	tv_start = usage.ru_utime;
#endif
}


/**
 * Measures process time (CPU time) elapsed between the last call to start()
 * and the time this method is called.
 *
 * @return elapsed process time in seconds.
 */
inline double ProcessStopwatch::elapsed() const {
#if defined(__unix) && defined(__sun)
	hrtime_t now = gethrvtime();
	// convert to seconds since gethrvtime() returns nanoseconds
	return (now-start_time)/1e9;
#elif defined(__linux__)
	rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	timeval tv_end = usage.ru_utime;
	timeval tv_elapsed;
	timersub(&tv_end, &tv_start, &tv_elapsed);
	return double(tv_elapsed.tv_sec) + tv_elapsed.tv_usec * 1e-6;
#endif
}
