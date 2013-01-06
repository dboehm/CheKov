package logger;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;

public class LogReadConfigExample {

	private static final long FIVE_SECS = 5_000;
	private static final Logger logger = Logger.getLogger(LogReadConfigExample.class);
	private static final long TWO_SECS = 2_000;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		PropertyConfigurator.configureAndWatch("config/log4j.properties",
				FIVE_SECS);
		logger.info("LogReadConfigExample started");

		while (!Thread.currentThread().isInterrupted()) {
			logger.debug("DEBUG");
			logger.info("INFO");
			logger.warn("WARN");
			logger.error("ERROR");
			logger.fatal("FATAL");
			try {
				Thread.sleep(TWO_SECS);
			} catch (final InterruptedException e) {
				Thread.currentThread().interrupt();
			}
		}

	}

}
