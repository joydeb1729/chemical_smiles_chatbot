"""
Utility functions for logging.
"""

import os
import logging
from logging.handlers import RotatingFileHandler

from config.settings import LOG_LEVEL, LOG_FILE

def setup_logger(name, level=None, log_file=None):
    """
    Set up a logger with the specified name.
    
    Args:
        name (str): The name of the logger.
        level (str, optional): The logging level. Defaults to config value.
        log_file (str, optional): The log file path. Defaults to config value.
        
    Returns:
        logging.Logger: The configured logger.
    """
    # Get the logger
    logger = logging.getLogger(name)
    
    # Skip if the logger is already configured
    if logger.handlers:
        return logger
    
    # Set log level
    level = level or LOG_LEVEL
    log_level = getattr(logging, level.upper())
    logger.setLevel(log_level)
    
    # Create formatters
    console_formatter = logging.Formatter('%(levelname)s - %(message)s')
    file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    # Create file handler if log file is specified
    log_file = log_file or LOG_FILE
    if log_file:
        # Create log directory if it doesn't exist
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir, exist_ok=True)
        
        # Create rotating file handler
        file_handler = RotatingFileHandler(
            log_file,
            maxBytes=10 * 1024 * 1024,  # 10 MB
            backupCount=5
        )
        file_handler.setLevel(log_level)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
    
    return logger
