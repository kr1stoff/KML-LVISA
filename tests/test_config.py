from kml_lvisa import get_threads_dict, get_conda_env_dict, get_database_dict, get_my_scripts_path, get_software_dict


def test_threads():
    print(get_threads_dict())


def test_conda_env():
    print(get_conda_env_dict())


def test_database():
    print(get_database_dict())


def test_my_scripts_path():
    print(get_my_scripts_path())


def test_software():
    print(get_software_dict())
