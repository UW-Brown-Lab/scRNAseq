  // Function to toggle visibility of elements with class "hljs"
  function toggleHighlightElements() {
    var elements = document.querySelectorAll('.fold-hide');
    elements.forEach(function(element) {
      if (element.style.display === 'none') {
        element.style.display = 'block';
      } else {
        element.style.display = 'none';
      }
    });
  }
  // Function to create and append a button for toggling visibility
  function createToggleButton() {
    var button = document.createElement('button');
    button.textContent = 'Toggle Code Visibility';
    button.addEventListener('click', toggleHighlightElements);
    button.style.backgroundColor = 'white';
    button.style.border = '2px solid black';
    button.style.color = 'black';
    button.style.padding = '15px 32px';
    button.style.textAlign = 'center';
    button.style.textDecoration = 'none';
    button.style.display = 'inline-block';
    button.style.fontSize = '16px';
    button.style.fontFamily = 'et-book, Palatino, "Palatino Linotype", "Palatino LT STD", "Book Antiqua", Georgia, serif';
    var header = document.querySelector('.date'); // Change selector to match your header element
    var buttonDiv = document.createElement('div');
    buttonDiv.id = 'button-div';
    header.insertAdjacentElement('afterend',buttonDiv);
    var divElement = document.getElementById('button-div');
    divElement.appendChild(button);
  }
  // Call the function to create and append the toggle button when the document is loaded
  document.addEventListener('DOMContentLoaded', createToggleButton);
  document.addEventListener('DOMContentLoaded', toggleHighlightElements);