from kml_lvisa import get_threads_dict, get_conda_env_dict, get_database_dict, get_my_scripts_path, get_software_dict
from pprint import pprint

print("Threads")
print(get_threads_dict())

print("Conda Environment")
print(get_conda_env_dict())

print("Database")
pprint(get_database_dict())

print("My Scripts Path")
print(get_my_scripts_path())

print("Software")
print(get_software_dict())
