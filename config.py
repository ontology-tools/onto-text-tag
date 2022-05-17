import os

# Example configuration
DEBUG = True

APP_TITLE = "Ontology Text Tagging"

# DATABASE_URI = 'sqlite:////tmp/github-flask-ontospreaded.db'

if os.environ.get("FLASK_ENV")=='development':
    SECRET_KEY = os.environ.get('FLASK_SECRET_KEY')
else:
     # Import the Secret Manager client library.
     #from google.cloud import secretmanager
     # Create the Secret Manager client.
     #client = secretmanager.SecretManagerServiceClient()

     #project_id = "onto-text-tag"
#     client_id = "GITHUB_CLIENT_ID"
#     # Build the resource name of the secret version.
#     name = f"projects/{project_id}/secrets/{client_id}/versions/latest"
#     # Access the secret version.
#     response = client.access_secret_version(request={"name": name})
#     # Build the resource name of the secret version.
#     name = f"projects/{project_id}/secrets/{client_secret}/versions/latest"
#     # Access the secret version.
#     response = client.access_secret_version(request={"name": name})

     #flask_secret = "FLASK_SECRET_KEY"
     # Build the resource name of the secret version.
     #name = f"projects/{project_id}/secrets/{flask_secret}/versions/latest"
     # Access the secret version.
     #response = client.access_secret_version(request={"name": name})
     #SECRET_KEY = response.payload.data.decode("UTF-8")
     # SECRET_KEY = os.environ.get('FLASK_SECRET_KEY')
     SECRET_KEY = "test_static_key" #todo: change to environment variable


