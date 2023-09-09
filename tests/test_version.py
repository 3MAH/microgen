import microgen

def test_version():
    package_version = microgen.__version__

    with open(file="conda.recipe/meta.yaml", mode="r", encoding='utf-8') as file:
        line = file.readline()
        conda_package_version = line.split('"')[1]

    assert package_version == conda_package_version
