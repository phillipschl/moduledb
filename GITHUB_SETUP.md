# Setting Up GitHub Repository

Follow these steps to link this project to your GitHub account:

## 1. Initialize Git Repository

If you haven't already initialized this as a Git repository, run the following commands in the project directory:

```bash
git init
git add .
git commit -m "Initial commit"
```

## 2. Create GitHub Repository

1. Go to [GitHub](https://github.com/) and sign in to your account
2. Click on the "+" icon in the top-right corner and select "New repository"
3. Name the repository "moduledb"
4. Optionally add a description: "A pipeline for generating a database of NRPS C-A-T modules"
5. Leave the repository as Public (or select Private if preferred)
6. Do NOT initialize the repository with a README, .gitignore, or license (as we already have these files locally)
7. Click "Create repository"

## 3. Link Local Repository to GitHub

After creating the repository, GitHub will show a page with setup instructions. Use the commands under "…or push an existing repository from the command line":

```bash
git remote add origin https://github.com/phillipschl/moduledb.git
git branch -M main
git push -u origin main
```

## 4. Verify Setup

1. Refresh your GitHub repository page
2. You should see all your project files now available on GitHub

## 5. Future Updates

For future changes to your code, use the standard Git workflow:

```bash
# Make changes to your files
git add .                         # Stage all changes
git commit -m "Description of changes"   # Commit changes with a message
git push                          # Push changes to GitHub
``` 