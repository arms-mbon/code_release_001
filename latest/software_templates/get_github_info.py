import requests
import json
import base64

# SIDENOTE: only public repositories
owner = "vliz-be-opsci"
repo = "k-gap"

urls = {"github": f"https://api.github.com/repos/{owner}/{repo}", 
        "contributors": f"https://api.github.com/repos/{owner}/{repo}/contributors",
        "languages": f"https://api.github.com/repos/{owner}/{repo}/languages"
        }

for name,url in urls.items():
    response = requests.get(url)
    if response.status_code == 200:
        print(response.json())
        with open(f'{owner}_{repo}_{name}.json', 'w') as json_file:
            json.dump(response.json(), json_file, indent=4)

def get_file_contents(owner, repo, path):
    url = f"https://api.github.com/repos/{owner}/{repo}/contents/{path}"
    response = requests.get(url)
    if response.status_code == 200:
        content = response.json()
        if 'content' in content:
            file_content = base64.b64decode(content['content']).decode('utf-8')
            return file_content
        else:
            return None
    else:
        return None


# List of potential requirements files
requirements_files = ['requirements.txt', 'package.json', 'Gemfile', 'Pipfile', 'docker-compose.yml']

requirements_data = {}

for file in requirements_files:
    content = get_file_contents(owner, repo, file)
    if content:
        requirements_data[file] = content

# Write the requirements data to a JSON file
if requirements_data:
    with open(f'{owner}_{repo}_requirements.json', 'w') as json_file:
        json.dump(requirements_data, json_file, indent=4)
    print("Requirements information has been written to 'requirements.json' file.")
else:
    print("No requirements files found in the repository.")