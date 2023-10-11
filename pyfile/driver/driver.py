import numpy as np
import func
import sys

driver = func.DRIVER()

if(len(sys.argv) > 1):
    if(sys.argv[1] == 'clean'):
        driver.delete_files()

    if(sys.argv[1] == 'create_rules'):
        driver.create_rules()

    if(sys.argv[1] == 'run'):
        driver.create_rules()
        if(len(sys.argv) > 4):
            driver.current_thread = int(sys.argv[2])
            driver.num_thread = int(sys.argv[3])
            driver.cuda_device = int(sys.argv[4])
        driver.run()

    if(sys.argv[1] == 'run_rules'):
        if(len(sys.argv) > 4):
            driver.current_thread = int(sys.argv[2])
            driver.num_thread = int(sys.argv[3])
            driver.cuda_device = int(sys.argv[4])
        driver.run()



























#