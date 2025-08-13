/**
 * JavaScript functions for the Chemical SMILES Chatbot application.
 */

document.addEventListener('DOMContentLoaded', function() {
    // Initialize the application
    initApp();
});

/**
 * Initialize the application.
 */
function initApp() {
    // Set up the form submit handler
    const queryForm = document.getElementById('query-form');
    if (queryForm) {
        queryForm.addEventListener('submit', handleQuerySubmit);
    }
    
    // Set up the clear button handler
    const clearButton = document.getElementById('clear-button');
    if (clearButton) {
        clearButton.addEventListener('click', clearResults);
    }
}

/**
 * Handle the query form submission.
 * 
 * @param {Event} event - The form submit event.
 */
function handleQuerySubmit(event) {
    // Prevent the default form submission
    event.preventDefault();
    
    // Get the query value
    const queryInput = document.getElementById('query-input');
    const query = queryInput.value.trim();
    
    if (!query) {
        showError('Please enter a chemical query.');
        return;
    }
    
    // Show loading spinner
    showLoading(true);
    
    // In a real application, this would submit the form via AJAX
    // Here we'll just simulate it
    setTimeout(function() {
        // Hide loading spinner
        showLoading(false);
        
        // Show success message
        showSuccess('Query processed successfully!');
    }, 1500);
}

/**
 * Show or hide the loading spinner.
 * 
 * @param {boolean} isLoading - Whether to show the loading spinner.
 */
function showLoading(isLoading) {
    const loadingSpinner = document.getElementById('loading-spinner');
    if (loadingSpinner) {
        loadingSpinner.style.display = isLoading ? 'flex' : 'none';
    }
}

/**
 * Show an error message.
 * 
 * @param {string} message - The error message to show.
 */
function showError(message) {
    const errorBox = document.getElementById('error-box');
    if (errorBox) {
        errorBox.textContent = message;
        errorBox.style.display = 'block';
        
        // Hide the error after 5 seconds
        setTimeout(function() {
            errorBox.style.display = 'none';
        }, 5000);
    }
}

/**
 * Show a success message.
 * 
 * @param {string} message - The success message to show.
 */
function showSuccess(message) {
    const successBox = document.getElementById('success-box');
    if (successBox) {
        successBox.textContent = message;
        successBox.style.display = 'block';
        
        // Hide the success message after 5 seconds
        setTimeout(function() {
            successBox.style.display = 'none';
        }, 5000);
    }
}

/**
 * Clear the results.
 */
function clearResults() {
    const resultsContainer = document.getElementById('results-container');
    if (resultsContainer) {
        resultsContainer.innerHTML = '';
    }
    
    const queryInput = document.getElementById('query-input');
    if (queryInput) {
        queryInput.value = '';
    }
}
