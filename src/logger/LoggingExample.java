package logger;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

public class LoggingExample {

	private static final Logger LOGGER = Logger.getRootLogger();
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// Default configuration to console
		BasicConfigurator.configure();
		
		// set Level
		LOGGER.setLevel(Level.ERROR);
		// Write Log messages
		LOGGER.info("Info message from LoggingExample");
		LOGGER.error("Error-message from Logging Example");
		

	}

}
