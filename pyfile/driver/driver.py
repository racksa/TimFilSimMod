import numpy as np
import func

driver = func.DRIVER('config.ini')
driver.run()

# commands = [
#         "ls -l",
#         "echo 'Hello'"
#     ]
# driver.run_commands(commands)