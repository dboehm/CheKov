package logger;

import java.io.File;
import java.io.IOException;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;

public class LoggingConfigExample {

	private static final Logger LOGGER = Logger
			.getLogger(LoggingConfigExample.class);
	private static final String LOGFILE = "LogFile.log";
	private static final boolean APPEND = true;

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		// Create layout
		final SimpleLayout layout = new SimpleLayout();

		// Create and attach ConsoleAppender to the SimpleLayout
		LOGGER.addAppender(new ConsoleAppender(layout));

		// Create and attach FileAppender to the logger

		try {
			LOGGER.addAppender(new FileAppender(layout, LOGFILE, APPEND));
		} catch (final IOException ex) {
			LOGGER.warn(
					"Can't create FileAppendar for "
							+ new File(LOGFILE).getAbsolutePath(), ex);
		}

		// set Level
		LOGGER.setLevel(Level.DEBUG);
		// Write Log messages
		LOGGER.info("Info message from LoggingExample"); // Filtered
		LOGGER.warn("Warn message from LoggingExample");
		LOGGER.error("Error-message from Logging Example");
		LOGGER.debug("Debugger says debugging");

	}
}
